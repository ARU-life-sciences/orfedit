//! Minimal multi-FASTA loading (whole per-species file, once, reused across
//! every gene scanned for that species - these are single mitogenome/
//! plastome assemblies, a few hundred kb to a couple Mb, trivially fine to
//! hold in memory) plus the genomic-window/strand/frame <-> codon-index
//! coordinate mapping the aligner's output needs translated back through.

use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

pub struct FastaIndex {
    pub sequences: HashMap<String, Vec<u8>>,
}

impl FastaIndex {
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        let mut sequences = HashMap::new();
        let mut cur_id: Option<String> = None;
        let mut cur_seq: Vec<u8> = Vec::new();
        for line in text.lines() {
            if let Some(id) = line.strip_prefix('>') {
                if let Some(prev_id) = cur_id.take() {
                    sequences.insert(prev_id, std::mem::take(&mut cur_seq));
                }
                cur_id = Some(id.split_whitespace().next().unwrap_or("").to_string());
            } else {
                cur_seq.extend(line.trim_end().bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        if let Some(prev_id) = cur_id.take() {
            sequences.insert(prev_id, cur_seq);
        }
        Ok(FastaIndex { sequences })
    }

    /// 0-based half-open [start, end) on the given contig, sense-corrected
    /// for strand ('+' or '-'). Clamped to contig bounds.
    pub fn extract_window(&self, contig: &str, start: i64, end: i64, strand: char) -> Result<(usize, usize, Vec<u8>)> {
        let seq = self
            .sequences
            .get(contig)
            .with_context(|| format!("contig {:?} not found in fasta", contig))?;
        let clamped_start = start.max(0) as usize;
        let clamped_end = (end.max(0) as usize).min(seq.len());
        if clamped_start >= clamped_end {
            bail!("empty window [{}, {}) on contig {} (len {})", start, end, contig, seq.len());
        }
        let mut slice = seq[clamped_start..clamped_end].to_vec();
        if strand == '-' {
            revcomp_in_place(&mut slice);
        }
        Ok((clamped_start, clamped_end, slice))
    }
}

pub fn revcomp_in_place(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = complement(*b);
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other, // N and anything else passes through
    }
}

/// Split a sense-oriented window sequence into complete codons for one of
/// the 3 reading frames (frame in {0,1,2} = how many leading bases to skip).
pub fn codons_for_frame(sense_seq: &[u8], frame: usize) -> Vec<[u8; 3]> {
    sense_seq[frame..]
        .chunks_exact(3)
        .map(|c| [c[0], c[1], c[2]])
        .collect()
}

/// Map a codon-index span [first, last_exclusive) in one frame's codon list
/// back to absolute genomic 0-based half-open [start, end) coordinates,
/// given the window this frame's codons were cut from.
///
/// `window_start`/`window_end` are the *genomic* (always-forward-strand)
/// bounds the sense sequence was extracted from (i.e. what
/// `FastaIndex::extract_window` returned as clamped_start/clamped_end,
/// before any revcomp).
pub fn codon_span_to_genomic(
    window_start: usize,
    window_end: usize,
    strand: char,
    frame: usize,
    first_codon: usize,
    last_codon_exclusive: usize,
) -> (usize, usize) {
    sense_span_to_genomic(window_start, window_end, strand, frame + first_codon * 3, frame + last_codon_exclusive * 3)
}

/// Map a sense-relative half-open span [sense_start, sense_end) of the
/// window back to absolute genomic 0-based half-open [start, end).
pub fn sense_span_to_genomic(
    window_start: usize,
    window_end: usize,
    strand: char,
    sense_start: usize,
    sense_end: usize,
) -> (usize, usize) {
    match strand {
        '+' => (window_start + sense_start, window_start + sense_end),
        // Sense position p <-> genomic index (window_end - 1 - p), so a
        // sense-forward span runs backwards on the genome.
        _ => (window_end - sense_end, window_end - sense_start),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plus_strand_roundtrip() {
        // window [100, 130), frame 0, codons 2..5 (sense bases 6..15)
        let (s, e) = codon_span_to_genomic(100, 130, '+', 0, 2, 5);
        assert_eq!((s, e), (106, 115));
    }

    #[test]
    fn minus_strand_roundtrip() {
        // A 12bp window, '-' strand, frame 0: sense codon 0 is the LAST 3
        // genomic bases; sense codon 3 (last, 4 codons total) is the FIRST
        // 3 genomic bases. Aligning sense codons [0,1) should map to the
        // genomic 3' end of the window on the + coordinate line.
        let (s, e) = codon_span_to_genomic(100, 112, '-', 0, 0, 1);
        assert_eq!((s, e), (109, 112));
        let (s2, e2) = codon_span_to_genomic(100, 112, '-', 0, 3, 4);
        assert_eq!((s2, e2), (100, 103));
    }

    #[test]
    fn revcomp_is_correct() {
        let mut seq = b"ACGTN".to_vec();
        revcomp_in_place(&mut seq);
        assert_eq!(seq, b"NACGT");
    }
}
