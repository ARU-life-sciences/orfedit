pub mod align;
pub mod codon;
pub mod fasta;
pub mod profile;

use align::NtColumnCall;
use anyhow::Result;
use codon::AA_ORDER;
use fasta::{sense_span_to_genomic, FastaIndex};
use profile::Profile;

pub struct EditCall {
    pub genomic_pos: usize,
    pub profile_column: usize, // 0-indexed
    pub resulting_aa: char,
    /// A U->C edit (the genome has T, the transcript reads C); false for
    /// the usual C->U.
    pub reverse: bool,
}

/// A frameshift in the corrected gene: genomic 0-based half-open
/// [start, end) of the 1 or 2 bases the alignment skipped to change frame,
/// lying between profile columns `profile_column - 1` and `profile_column`
/// (see `align::Frameshift`).
pub struct FrameshiftCall {
    pub start: usize,
    pub end: usize,
    pub n_bases: usize,
    pub profile_column: usize,
    /// Genomic [start, end) the frame change could equally well sit in
    /// (align::Frameshift::window); contains [start, end).
    pub window_start: usize,
    pub window_end: usize,
}

pub struct ScanOutcome {
    /// Reading frame (0-2, within the window) of the first aligned codon.
    pub frame: usize,
    pub score: f64,
    pub corrected_start: usize,
    pub corrected_end: usize,
    pub protein: String, // '-' for deleted (unaligned) profile columns
    pub edits: Vec<EditCall>,
    pub n_columns_aligned: usize,
    pub n_columns_total: usize,
    /// Whether a genuine (possibly edit-created) stop codon was found
    /// immediately after the aligned protein and folded into
    /// corrected_start/corrected_end. `profile` alignments are built from
    /// stop-stripped protein translations, so without this the corrected
    /// span would systematically land one codon short of `.ctg.bed`'s own
    /// stop-inclusive convention on essentially every call, regardless of
    /// whether real editing happened - confirmed empirically (see
    /// `analysis/editing/README.md` in the calling repo).
    pub includes_stop: bool,
    /// Frame changes inside the gene, in genomic order (empty unless
    /// `ScanParams::frameshift_penalty` is finite).
    pub frameshifts: Vec<FrameshiftCall>,
}

/// A single sense-codon position (bit 0/1/2 within a codon) mapped back to
/// its absolute genomic coordinate. Mirrors `codon_span_to_genomic`'s
/// range math for one base instead of a whole span.
fn sense_base_to_genomic(window_start: usize, window_end: usize, strand: char, sense_pos: usize) -> usize {
    match strand {
        '+' => window_start + sense_pos,
        _ => window_end - 1 - sense_pos,
    }
}

/// Optional filters on frame changes (both off by default). With
/// `frameshift_flank` > 0 a frame change is kept only when the alignment is
/// sound on both sides of it: that many profile columns each side, at most
/// `frameshift_max_deleted` of them unaligned, the aligned ones scoring at
/// least `frameshift_min_flank_score` on average. With `frameshift_min_gain`
/// > 0 a gene keeps its frame changes only if they beat its best
/// single-frame alignment by that much. If either fails the gene gets its
/// single-frame alignment. Measured on 1,500 1-bp indels planted in the
/// coding sequence of DToL plant organelle genes, the flank filter (10,
/// 2, 1.5) at penalty 30 cost 8 points of recall (90% -> 81%) for few
/// fewer one-off calls in good genes (57 -> 45 over 2,122 assemblies) -
/// most of what it removes is in partial copies or recurs across
/// assemblies - so neither is used by default.
pub const FRAMESHIFT_FLANK: usize = 0;
pub const FRAMESHIFT_MAX_DELETED: usize = 2;
pub const FRAMESHIFT_MIN_FLANK_SCORE: f64 = 1.5;
/// A gene keeps its frame changes only if its alignment with them beats
/// its best single-frame alignment by at least this much (after paying
/// the frameshift penalty).
pub const FRAMESHIFT_MIN_GAIN: f64 = 0.0;

/// Whether every frame change in `result` has sound alignment both sides
/// (always, with `params.frameshift_flank` 0).
impl ScanParams {
    /// The edit costs these parameters ask for.
    pub fn edit_costs(&self) -> codon::EditCosts {
        codon::EditCosts {
            c_to_u: self.edit_penalty,
            u_to_c: self.reverse_edit_penalty,
            u_to_c_stop: self.stop_edit_penalty,
            c_to_u_known: self.known_edit_penalty.unwrap_or(self.edit_penalty),
            c_to_u_novel: self.novel_edit_penalty.unwrap_or(self.edit_penalty),
        }
    }
}

pub fn frameshifts_supported(profile: &Profile, result: &align::NtAlignResult, params: &ScanParams) -> bool {
    let n = result.columns.len();
    let k = params.frameshift_flank;
    if k == 0 {
        return true;
    }
    result.frameshifts.iter().all(|f| {
        let c = f.after_column;
        if c < k || c + k > n {
            return false;
        }
        [c - k..c, c..c + k].into_iter().all(|range| {
            let (mut sum, mut aligned, mut deleted) = (0.0, 0usize, 0usize);
            for j in range {
                match &result.columns[j] {
                    NtColumnCall::Aligned { aa, .. } => {
                        sum += profile.score(j, *aa);
                        aligned += 1;
                    }
                    NtColumnCall::Deleted => deleted += 1,
                }
            }
            deleted <= params.frameshift_max_deleted
                && aligned > 0
                && sum / aligned as f64 >= params.frameshift_min_flank_score
        })
    })
}

pub struct ScanParams {
    pub gap_open: f64,
    pub gap_extend: f64,
    pub edit_penalty: f64,
    pub flank: i64,
    /// Cost of skipping 1 or 2 bases to change reading frame inside a gene;
    /// `f64::INFINITY` = single-frame alignment only.
    pub frameshift_penalty: f64,
    /// Profile columns each side of a frame change that must align soundly
    /// (0: keep every frame change), how many of them may be unaligned, and
    /// their minimum mean score - see FRAMESHIFT_FLANK.
    pub frameshift_flank: usize,
    pub frameshift_max_deleted: usize,
    pub frameshift_min_flank_score: f64,
    /// See FRAMESHIFT_MIN_GAIN.
    pub frameshift_min_gain: f64,
    /// When no stop codon directly follows the aligned protein, look this
    /// many codons further on, in frame, for a literal one (TAA/TAG/TGA)
    /// and extend the gene to it (0: don't). Only when the alignment
    /// reaches the profile's end, give or take this many columns. A
    /// profile often lacks a gene's last residue or two (boundary-truncated
    /// training calls), which left the corrected gene those codons and the
    /// stop short.
    pub stop_search: usize,
    /// Cost of a U->C ("reverse") edit that is not a stop read-through;
    /// infinite (the default) for none. Ferns, some lycophytes and
    /// hornworts edit both ways.
    pub reverse_edit_penalty: f64,
    /// Cost of reading a genomic in-frame stop codon as a residue by a
    /// U->C edit (TAA, TAG -> Gln; TGA -> Arg); infinite (the default) for
    /// none. Kept apart from `reverse_edit_penalty` because a stop inside
    /// a gene is far stronger evidence of an edit than any other T.
    pub stop_edit_penalty: f64,
    /// For a profile that says where its reference proteins are edited
    /// (its `edit` column): the cost of a C->U edit at a column edited in
    /// them, and at one that is not. `None` (the default) for
    /// `edit_penalty`, which is also the cost at every column of a profile
    /// without that column.
    pub known_edit_penalty: Option<f64>,
    pub novel_edit_penalty: Option<f64>,
}

impl Default for ScanParams {
    fn default() -> Self {
        ScanParams {
            gap_open: 8.0,
            gap_extend: 1.0,
            edit_penalty: 2.5,
            flank: 300,
            frameshift_penalty: f64::INFINITY,
            frameshift_flank: FRAMESHIFT_FLANK,
            frameshift_max_deleted: FRAMESHIFT_MAX_DELETED,
            frameshift_min_flank_score: FRAMESHIFT_MIN_FLANK_SCORE,
            frameshift_min_gain: FRAMESHIFT_MIN_GAIN,
            stop_search: 0,
            reverse_edit_penalty: f64::INFINITY,
            stop_edit_penalty: f64::INFINITY,
            known_edit_penalty: None,
            novel_edit_penalty: None,
        }
    }
}

/// Scan one (species, gene) anchor: extract a flanked genomic window around
/// `anchor_start`/`anchor_end` (0-based half-open, the existing `.ctg.bed`
/// call), try all 3 reading frames against `profile`, and return the
/// best-scoring edit-tolerant alignment.
pub fn scan_gene(
    fasta_idx: &FastaIndex,
    contig: &str,
    anchor_start: i64,
    anchor_end: i64,
    strand: char,
    profile: &Profile,
    params: &ScanParams,
) -> Result<ScanOutcome> {
    let (window_start, window_end, sense_seq) = fasta_idx.extract_window(
        contig,
        anchor_start - params.flank,
        anchor_end + params.flank,
        strand,
    )?;

    let mut result = align::align_nt(
        profile,
        &sense_seq,
        params.gap_open,
        params.gap_extend,
        params.edit_costs(),
        params.frameshift_penalty,
    );
    let filtering = params.frameshift_min_gain > 0.0 || params.frameshift_flank > 0;
    if filtering && !result.frameshifts.is_empty() {
        let single = align::align_nt(profile, &sense_seq, params.gap_open, params.gap_extend, params.edit_costs(), f64::INFINITY);
        if result.score - single.score < params.frameshift_min_gain || !frameshifts_supported(profile, &result, params) {
            result = single;
        }
    }
    let frame = result.sense_span.map_or(0, |(s, _)| s % 3);

    let (mut corrected_start, mut corrected_end) = match result.sense_span {
        Some((s, e)) => sense_span_to_genomic(window_start, window_end, strand, s, e),
        None => (anchor_start.max(0) as usize, anchor_end.max(0) as usize),
    };

    let mut protein = String::with_capacity(result.columns.len());
    let mut edits = Vec::new();
    let mut n_aligned = 0usize;
    let mut includes_stop = false;
    for (col_idx, call) in result.columns.iter().enumerate() {
        match call {
            NtColumnCall::Aligned { sense_start, aa, n_edits, edit_mask } => {
                n_aligned += 1;
                protein.push(*aa as char);
                if *n_edits > 0 {
                    for bit in 0..3u8 {
                        if edit_mask & (1 << bit) != 0 {
                            let sense_pos = sense_start + bit as usize;
                            let genomic_pos = sense_base_to_genomic(window_start, window_end, strand, sense_pos);
                            // the genome's base says which way: a C was read as U, a T as C
                            let reverse = sense_seq[sense_pos] == b'T';
                            edits.push(EditCall { genomic_pos, profile_column: col_idx, resulting_aa: *aa as char, reverse });
                        }
                    }
                }
            }
            NtColumnCall::Deleted => protein.push('-'),
        }
    }
    // Fold in a trailing stop codon, if one is genuinely there (possibly
    // only after a C->T edit) - see `ScanOutcome::includes_stop`.
    if let Some((_, sense_end)) = result.sense_span {
        if sense_end + 3 <= sense_seq.len() {
            let stop_codon = [sense_seq[sense_end], sense_seq[sense_end + 1], sense_seq[sense_end + 2]];
            if let Some(hyp) = codon::codon_hypotheses(stop_codon)
                .into_iter()
                .filter(|h| h.aa == codon::STOP)
                .min_by_key(|h| h.n_edits)
            {
                // A stop that needs an edit (CAA/CAG/CGA), with a literal
                // one straight after it: the codon is the protein's last
                // residue (Gln/Arg, missing from the profile) and the
                // literal stop the gene's.
                let p = sense_end + 3;
                let residue_then_stop = hyp.n_edits > 0
                    && params.stop_search > 0
                    && p + 3 <= sense_seq.len()
                    && codon::translate_codon([sense_seq[p], sense_seq[p + 1], sense_seq[p + 2]]) == codon::STOP;
                let stop_end = if residue_then_stop { p + 3 } else { p };
                let (ext_start, ext_end) = sense_span_to_genomic(window_start, window_end, strand, sense_end, stop_end);
                if strand == '+' {
                    corrected_end = ext_end;
                } else {
                    corrected_start = ext_start;
                }
                includes_stop = true;
                if hyp.n_edits > 0 && !residue_then_stop {
                    for bit in 0..3u8 {
                        if hyp.edit_mask & (1 << bit) != 0 {
                            let sense_pos = sense_end + bit as usize;
                            let genomic_pos = sense_base_to_genomic(window_start, window_end, strand, sense_pos);
                            edits.push(EditCall { genomic_pos, profile_column: result.columns.len(), resulting_aa: '*', reverse: false });
                        }
                    }
                }
            }
        }
    }
    // No stop right there: a literal one a few codons on, with the
    // alignment at (or nearly at) the profile's end, is the gene's own.
    if !includes_stop && params.stop_search > 0 {
        let trailing_deleted =
            result.columns.iter().rev().take_while(|c| matches!(c, NtColumnCall::Deleted)).count();
        if let (Some((_, sense_end)), true) = (result.sense_span, trailing_deleted <= params.stop_search) {
            let stop_at = (1..=params.stop_search).map(|i| sense_end + 3 * i).find(|&p| {
                p + 3 <= sense_seq.len()
                    && codon::translate_codon([sense_seq[p], sense_seq[p + 1], sense_seq[p + 2]]) == codon::STOP
            });
            if let Some(p) = stop_at {
                let (ext_start, ext_end) = sense_span_to_genomic(window_start, window_end, strand, sense_end, p + 3);
                if strand == '+' {
                    corrected_end = ext_end;
                } else {
                    corrected_start = ext_start;
                }
                includes_stop = true;
            }
        }
    }
    let mut frameshifts: Vec<FrameshiftCall> = result
        .frameshifts
        .iter()
        .map(|f| {
            let (start, end) =
                sense_span_to_genomic(window_start, window_end, strand, f.sense_pos, f.sense_pos + f.n_bases);
            let (ws, we) = sense_span_to_genomic(window_start, window_end, strand, f.window.0, f.window.1);
            FrameshiftCall { start, end, n_bases: f.n_bases, profile_column: f.after_column, window_start: ws, window_end: we }
        })
        .collect();
    frameshifts.sort_by_key(|f| f.start);
    edits.sort_by_key(|e| e.genomic_pos);

    Ok(ScanOutcome {
        frame,
        score: result.score,
        corrected_start,
        corrected_end,
        protein,
        edits,
        n_columns_aligned: n_aligned,
        n_columns_total: profile.len(),
        includes_stop,
        frameshifts,
    })
}

pub fn amino_acid_alphabet() -> &'static [u8; 20] {
    &AA_ORDER
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn mkv_profile() -> Profile {
        let mut columns = Vec::new();
        for target in [b'M', b'K', b'V'] {
            let mut col = [-5.0f64; 20];
            let idx = AA_ORDER.iter().position(|&x| x == target).unwrap();
            col[idx] = 5.0;
            columns.push(col);
        }
        Profile { gene: "test".to_string(), columns, edited: None }
    }

    fn fasta_with(contig: &str, seq: &str) -> FastaIndex {
        let mut sequences = HashMap::new();
        sequences.insert(contig.to_string(), seq.bytes().collect());
        FastaIndex { sequences }
    }

    #[test]
    fn literal_stop_is_folded_into_corrected_end() {
        let idx = fasta_with("ctg", "ATGAAAGTTTAA");
        let profile = mkv_profile();
        let params = ScanParams { flank: 3, ..ScanParams::default() };
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &params).unwrap();
        assert!(out.includes_stop);
        assert_eq!(out.corrected_end, 12);
        assert!(out.edits.is_empty());
    }

    #[test]
    fn edit_created_stop_is_folded_in_and_recorded() {
        // CAA (Gln, unedited) -> TAA (stop) via a single C->T edit.
        let idx = fasta_with("ctg", "ATGAAAGTTCAA");
        let profile = mkv_profile();
        let params = ScanParams { flank: 3, ..ScanParams::default() };
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &params).unwrap();
        assert!(out.includes_stop);
        assert_eq!(out.corrected_end, 12);
        assert_eq!(out.edits.len(), 1);
        assert_eq!(out.edits[0].resulting_aa, '*');
        assert_eq!(out.edits[0].genomic_pos, 9); // the C in CAA
    }

    #[test]
    fn non_stop_trailing_codon_is_not_folded_in() {
        let idx = fasta_with("ctg", "ATGAAAGTTGGG"); // GGG = Gly, never a stop
        let profile = mkv_profile();
        let params = ScanParams { flank: 3, ..ScanParams::default() };
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &params).unwrap();
        assert!(!out.includes_stop);
        assert_eq!(out.corrected_end, 9);
    }

    #[test]
    fn a_stop_a_few_codons_on_is_found_with_stop_search() {
        // profile M K V; the gene has two more residues (GGG CCC) before TAA
        let idx = fasta_with("ctg", "ATGAAAGTTGGGCCCTAAGGG");
        let profile = mkv_profile();
        let params = ScanParams { flank: 12, ..ScanParams::default() };
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &params).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end), (false, 9));
        let params = ScanParams { stop_search: 2, ..params };
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &params).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end), (true, 18));
        assert_eq!(out.protein, "MKV");
        // too far for a one-codon search
        let params = ScanParams { stop_search: 1, ..params };
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &params).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end), (false, 9));
        // CAA could be an edit-made stop, but TAA follows it: Gln, then the stop
        let idx = fasta_with("ctg", "ATGAAAGTTCAATAAGGG");
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &params).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end, out.edits.len()), (true, 15, 0));
        let out = scan_gene(&idx, "ctg", 0, 9, '+', &profile, &ScanParams { stop_search: 0, ..params }).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end, out.edits.len()), (true, 12, 1));
        // minus strand
        let idx = fasta_with("ctg", "CCCTTAGGGCCCAACTTTCAT");
        let params = ScanParams { stop_search: 2, ..params };
        let out = scan_gene(&idx, "ctg", 12, 21, '-', &profile, &params).unwrap();
        assert_eq!((out.includes_stop, out.corrected_start, out.corrected_end), (true, 3, 21));
    }

    fn one_residue_profile(residues: &[u8]) -> Profile {
        let mut columns = Vec::new();
        for &target in residues {
            let mut col = [-5.0f64; 20];
            col[AA_ORDER.iter().position(|&x| x == target).unwrap()] = 5.0;
            columns.push(col);
        }
        Profile { gene: "test".to_string(), columns, edited: None }
    }

    #[test]
    fn an_edit_at_a_column_the_references_do_not_edit_costs_more() {
        // profile M L K L F; the genome has Ser (TCA) at both Leu columns
        let mut profile = one_residue_profile(b"MLKLF");
        let idx = fasta_with("ctg", "GGGATGTCAAAATCATTTTAAGGG");
        let params = ScanParams { flank: 3, ..ScanParams::default() };
        let out = scan_gene(&idx, "ctg", 3, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.protein.as_str(), out.edits.len()), ("MLKLF", 2));
        // the profile says its references are edited at the first Leu only: without
        // separate costs that changes nothing
        profile.edited = Some(vec![false, true, false, false, false]);
        let out = scan_gene(&idx, "ctg", 3, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.protein.as_str(), out.edits.len()), ("MLKLF", 2));
        // with one for the other columns, the edit at the other Leu is no longer worth its cost
        let params = ScanParams { novel_edit_penalty: Some(20.0), ..params };
        let out = scan_gene(&idx, "ctg", 3, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.protein.as_str(), out.edits.len()), ("MLKSF", 1));
        assert_eq!((out.edits[0].genomic_pos, out.edits[0].profile_column, out.edits[0].resulting_aa), (7, 1, 'L'));
        // an edit too dear at the plain cost is still made at a known column, at its own cost
        let params = ScanParams { edit_penalty: 20.0, known_edit_penalty: Some(1.0), novel_edit_penalty: None, ..params };
        let out = scan_gene(&idx, "ctg", 3, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.protein.as_str(), out.edits.len()), ("MLKSF", 1));
        // and a profile that does not say takes the plain cost everywhere
        profile.edited = None;
        let out = scan_gene(&idx, "ctg", 3, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.protein.as_str(), out.edits.len()), ("MSKSF", 0));
    }

    #[test]
    fn reverse_editing_reads_an_in_frame_stop_as_a_residue() {
        // profile M K Q V F: the genome has TAA where the protein has Gln
        let profile = one_residue_profile(b"MKQVF");
        let idx = fasta_with("ctg", "GGGATGAAATAAGTTTTTTAGGGG");
        let params = ScanParams { flank: 3, ..ScanParams::default() };
        // C->U only: the stop can't be a residue, so its column is skipped
        let out = scan_gene(&idx, "ctg", 3, 18, '+', &profile, &params).unwrap();
        assert_ne!(out.protein, "MKQVF");
        assert!(out.edits.iter().all(|e| !e.reverse));
        // stop read-through alone is enough: TAA -> CAA, one U->C edit at the T
        let params = ScanParams { stop_edit_penalty: 2.5, ..params };
        let out = scan_gene(&idx, "ctg", 3, 18, '+', &profile, &params).unwrap();
        assert_eq!(out.protein, "MKQVF");
        assert_eq!(out.edits.len(), 1);
        assert_eq!((out.edits[0].genomic_pos, out.edits[0].resulting_aa, out.edits[0].reverse), (9, 'Q', true));
        // the literal stop after the last residue is still the gene's end
        assert!(out.includes_stop);
        assert_eq!((out.corrected_start, out.corrected_end), (3, 21));
        // minus strand: the same gene reverse-complemented
        let idx = fasta_with("ctg", "CCCCTAAAAAACTTATTTCATCCC");
        let out = scan_gene(&idx, "ctg", 6, 21, '-', &profile, &params).unwrap();
        assert_eq!(out.protein, "MKQVF");
        assert_eq!((out.edits.len(), out.edits[0].genomic_pos, out.edits[0].reverse), (1, 14, true));
        assert_eq!((out.corrected_start, out.corrected_end), (3, 21));
    }

    #[test]
    fn frameshift_is_reported_in_genomic_coordinates_on_both_strands() {
        // (M K V W F) x 6 with an extra A after the 15th codon, then a stop
        let half = "ATGAAAGTTTGGTTT".repeat(3);
        let fwd = format!("{half}A{half}TAA");
        let mut columns = Vec::new();
        for target in b"MKVWF".repeat(6).iter() {
            let mut col = [-5.0f64; 20];
            col[AA_ORDER.iter().position(|x| x == target).unwrap()] = 5.0;
            columns.push(col);
        }
        let profile = Profile { gene: "t".into(), columns, edited: None };
        let params = ScanParams { flank: 0, frameshift_penalty: 10.0, ..ScanParams::default() };
        let out = scan_gene(&fasta_with("c", &fwd), "c", 0, fwd.len() as i64, '+', &profile, &params).unwrap();
        assert_eq!(out.frameshifts.len(), 1);
        let f = &out.frameshifts[0];
        assert_eq!((f.start, f.end, f.n_bases, f.profile_column), (45, 46, 1, 15));
        assert_eq!((out.corrected_start, out.corrected_end, out.includes_stop), (0, 94, true));
        assert_eq!(out.protein, "MKVWF".repeat(6));

        // the same gene on the minus strand of a reverse-complemented contig
        let rev: String = fwd
            .bytes()
            .rev()
            .map(|b| match b { b'A' => 'T', b'T' => 'A', b'C' => 'G', _ => 'C' })
            .collect();
        let out = scan_gene(&fasta_with("c", &rev), "c", 0, rev.len() as i64, '-', &profile, &params).unwrap();
        let f = &out.frameshifts[0];
        assert_eq!((f.start, f.end, f.n_bases), (94 - 46, 94 - 45, 1));
        assert_eq!((out.corrected_start, out.corrected_end), (0, 94));

        // with the flank filter on, too near the gene's end to judge:
        // single-frame result instead
        let params = ScanParams { frameshift_flank: 10, frameshift_max_deleted: 2, ..params };
        let short = "ATGAAAGTTTGGTTTAATGAAAGTTTGGTTTTAA";
        let mut cols = Vec::new();
        for target in b"MKVWFMKVWF" {
            let mut col = [-5.0f64; 20];
            col[AA_ORDER.iter().position(|x| x == target).unwrap()] = 5.0;
            cols.push(col);
        }
        let p10 = Profile { gene: "t".into(), columns: cols, edited: None };
        let out = scan_gene(&fasta_with("c", short), "c", 0, short.len() as i64, '+', &p10, &params).unwrap();
        assert!(out.frameshifts.is_empty());
    }
}
