pub mod align;
pub mod codon;
pub mod fasta;
pub mod profile;

use align::ColumnCall;
use anyhow::Result;
use codon::AA_ORDER;
use fasta::{codon_span_to_genomic, codons_for_frame, FastaIndex};
use profile::Profile;

pub struct EditCall {
    pub genomic_pos: usize,
    pub profile_column: usize, // 0-indexed
    pub resulting_aa: char,
}

pub struct ScanOutcome {
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

pub struct ScanParams {
    pub gap_open: f64,
    pub gap_extend: f64,
    pub edit_penalty: f64,
    pub flank: i64,
}

impl Default for ScanParams {
    fn default() -> Self {
        ScanParams { gap_open: 8.0, gap_extend: 1.0, edit_penalty: 2.5, flank: 300 }
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

    let mut best: Option<(usize, align::AlignResult)> = None;
    for frame in 0..3 {
        let codons = codons_for_frame(&sense_seq, frame);
        let result = align::align(profile, &codons, params.gap_open, params.gap_extend, params.edit_penalty);
        let better = match &best {
            None => true,
            Some((_, prev)) => result.score > prev.score,
        };
        if better {
            best = Some((frame, result));
        }
    }
    let (frame, result) = best.expect("tried at least one frame");

    let (mut corrected_start, mut corrected_end) = match result.codon_span {
        Some((first, last_exclusive)) => {
            codon_span_to_genomic(window_start, window_end, strand, frame, first, last_exclusive)
        }
        None => (anchor_start.max(0) as usize, anchor_end.max(0) as usize),
    };

    let mut protein = String::with_capacity(result.columns.len());
    let mut edits = Vec::new();
    let mut n_aligned = 0usize;
    let mut includes_stop = false;
    for (col_idx, call) in result.columns.iter().enumerate() {
        match call {
            ColumnCall::Aligned { codon_idx, aa, n_edits, edit_mask } => {
                n_aligned += 1;
                protein.push(*aa as char);
                if *n_edits > 0 {
                    for bit in 0..3u8 {
                        if edit_mask & (1 << bit) != 0 {
                            let sense_pos = frame + codon_idx * 3 + bit as usize;
                            let genomic_pos = sense_base_to_genomic(window_start, window_end, strand, sense_pos);
                            edits.push(EditCall { genomic_pos, profile_column: col_idx, resulting_aa: *aa as char });
                        }
                    }
                }
            }
            ColumnCall::Deleted => protein.push('-'),
        }
    }
    // Fold in a trailing stop codon, if one is genuinely there (possibly
    // only after a C->T edit) - see `ScanOutcome::includes_stop`.
    if let Some((_, last_exclusive)) = result.codon_span {
        let next_sense_start = frame + last_exclusive * 3;
        if next_sense_start + 3 <= sense_seq.len() {
            let stop_codon = [
                sense_seq[next_sense_start],
                sense_seq[next_sense_start + 1],
                sense_seq[next_sense_start + 2],
            ];
            if let Some(hyp) = codon::codon_hypotheses(stop_codon)
                .into_iter()
                .filter(|h| h.aa == codon::STOP)
                .min_by_key(|h| h.n_edits)
            {
                let (ext_start, ext_end) =
                    codon_span_to_genomic(window_start, window_end, strand, frame, last_exclusive, last_exclusive + 1);
                if strand == '+' {
                    corrected_end = ext_end;
                } else {
                    corrected_start = ext_start;
                }
                includes_stop = true;
                if hyp.n_edits > 0 {
                    for bit in 0..3u8 {
                        if hyp.edit_mask & (1 << bit) != 0 {
                            let sense_pos = next_sense_start + bit as usize;
                            let genomic_pos = sense_base_to_genomic(window_start, window_end, strand, sense_pos);
                            edits.push(EditCall { genomic_pos, profile_column: result.columns.len(), resulting_aa: '*' });
                        }
                    }
                }
            }
        }
    }
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
}
