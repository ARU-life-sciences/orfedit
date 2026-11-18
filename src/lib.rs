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

/// A frame change is kept only when the alignment is sound on both sides
/// of it: FRAMESHIFT_FLANK profile columns each side, of which at most
/// FRAMESHIFT_MAX_DELETED have no codon, and the aligned ones score at
/// least FRAMESHIFT_MIN_FLANK_SCORE on average (profile score of the
/// amino acid read). Without this, a truncated or degenerate copy - where
/// part of the profile has nothing real to align to - picks up frame
/// changes that only let it score that junk in the best of three frames.
/// On the DToL plant organelle dataset (penalty 30) this cuts the genes
/// with frame changes from 1,973 to 698 while still finding 83% of 1-bp
/// indels planted in real genes (94% plastid, 72% mito). If any frame
/// change in a gene fails, the gene gets its single-frame alignment.
pub const FRAMESHIFT_FLANK: usize = 10;
pub const FRAMESHIFT_MAX_DELETED: usize = 2;
pub const FRAMESHIFT_MIN_FLANK_SCORE: f64 = 1.5;

/// Whether every frame change in `result` has sound alignment both sides.
fn frameshifts_supported(profile: &Profile, result: &align::NtAlignResult) -> bool {
    let n = result.columns.len();
    result.frameshifts.iter().all(|f| {
        let c = f.after_column;
        if c < FRAMESHIFT_FLANK || c + FRAMESHIFT_FLANK > n {
            return false;
        }
        [c - FRAMESHIFT_FLANK..c, c..c + FRAMESHIFT_FLANK].into_iter().all(|range| {
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
            deleted <= FRAMESHIFT_MAX_DELETED && aligned > 0 && sum / aligned as f64 >= FRAMESHIFT_MIN_FLANK_SCORE
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
}

impl Default for ScanParams {
    fn default() -> Self {
        ScanParams { gap_open: 8.0, gap_extend: 1.0, edit_penalty: 2.5, flank: 300, frameshift_penalty: f64::INFINITY }
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
        params.edit_penalty,
        params.frameshift_penalty,
    );
    if !result.frameshifts.is_empty() && !frameshifts_supported(profile, &result) {
        result = align::align_nt(profile, &sense_seq, params.gap_open, params.gap_extend, params.edit_penalty, f64::INFINITY);
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
                            edits.push(EditCall { genomic_pos, profile_column: col_idx, resulting_aa: *aa as char });
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
                let (ext_start, ext_end) =
                    sense_span_to_genomic(window_start, window_end, strand, sense_end, sense_end + 3);
                if strand == '+' {
                    corrected_end = ext_end;
                } else {
                    corrected_start = ext_start;
                }
                includes_stop = true;
                if hyp.n_edits > 0 {
                    for bit in 0..3u8 {
                        if hyp.edit_mask & (1 << bit) != 0 {
                            let sense_pos = sense_end + bit as usize;
                            let genomic_pos = sense_base_to_genomic(window_start, window_end, strand, sense_pos);
                            edits.push(EditCall { genomic_pos, profile_column: result.columns.len(), resulting_aa: '*' });
                        }
                    }
                }
            }
        }
    }
    let mut frameshifts: Vec<FrameshiftCall> = result
        .frameshifts
        .iter()
        .map(|f| {
            let (start, end) =
                sense_span_to_genomic(window_start, window_end, strand, f.sense_pos, f.sense_pos + f.n_bases);
            FrameshiftCall { start, end, n_bases: f.n_bases, profile_column: f.after_column }
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
        Profile { gene: "test".to_string(), columns }
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
        let profile = Profile { gene: "t".into(), columns };
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

        // too near the gene's end to judge: single-frame result instead
        let short = "ATGAAAGTTTGGTTTAATGAAAGTTTGGTTTTAA";
        let mut cols = Vec::new();
        for target in b"MKVWFMKVWF" {
            let mut col = [-5.0f64; 20];
            col[AA_ORDER.iter().position(|x| x == target).unwrap()] = 5.0;
            cols.push(col);
        }
        let p10 = Profile { gene: "t".into(), columns: cols };
        let out = scan_gene(&fasta_with("c", short), "c", 0, short.len() as i64, '+', &p10, &params).unwrap();
        assert!(out.frameshifts.is_empty());
    }
}
