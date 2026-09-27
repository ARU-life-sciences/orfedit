//! Standard genetic code (NCBI translation table 1) plus C->U (C->T at the
//! DNA level) edit-hypothesis enumeration.
//!
//! Plant mitochondrial/plastid RNA editing is (almost exclusively) C->U at
//! the transcript level, i.e. a genomic C can read as a T once transcribed
//! and edited. We never hypothesize other base changes.

/// Fixed amino acid column order used by `Profile` (`profile::AA_ORDER`
/// mirrors this exactly - keep them in sync).
pub const AA_ORDER: [u8; 20] = *b"ARNDCQEGHILKMFPSTWYV";
pub const STOP: u8 = b'*';
pub const UNKNOWN: u8 = b'X';

pub fn aa_index(aa: u8) -> Option<usize> {
    AA_ORDER.iter().position(|&x| x == aa)
}

/// Translate one codon (uppercase ACGT only; anything else -> UNKNOWN).
pub fn translate_codon(c: [u8; 3]) -> u8 {
    match c {
        [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => b'F',
        [b'T', b'T', b'A'] | [b'T', b'T', b'G'] => b'L',
        [b'T', b'C', _] => b'S',
        [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => b'Y',
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] => STOP,
        [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => b'C',
        [b'T', b'G', b'A'] => STOP,
        [b'T', b'G', b'G'] => b'W',
        [b'C', b'T', _] => b'L',
        [b'C', b'C', _] => b'P',
        [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => b'H',
        [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => b'Q',
        [b'C', b'G', _] => b'R',
        [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => b'I',
        [b'A', b'T', b'G'] => b'M',
        [b'A', b'C', _] => b'T',
        [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => b'N',
        [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => b'K',
        [b'A', b'G', b'T'] | [b'A', b'G', b'C'] => b'S',
        [b'A', b'G', b'A'] | [b'A', b'G', b'G'] => b'R',
        [b'G', b'T', _] => b'V',
        [b'G', b'C', _] => b'A',
        [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => b'D',
        [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => b'E',
        [b'G', b'G', _] => b'G',
        _ => UNKNOWN,
    }
}

/// One candidate reading of a genomic codon: the amino acid it produces and
/// how many C->T edits that reading assumes.
#[derive(Clone, Copy, Debug)]
pub struct CodonHypothesis {
    pub aa: u8,
    pub n_edits: u8,
    /// Bitmask (bits 0,1,2 = codon positions 0,1,2) of which positions were
    /// hypothesized edited.
    pub edit_mask: u8,
}

/// Every C->T reading of `codon` (2^k for k = number of Cs, k in 0..=3),
/// including the unedited (k=0) reading. Non-ACGT bases are left as-is (no
/// edit hypothesized there).
pub fn codon_hypotheses(codon: [u8; 3]) -> Vec<CodonHypothesis> {
    let c_positions: Vec<usize> = (0..3).filter(|&p| codon[p] == b'C').collect();
    let k = c_positions.len();
    let mut out = Vec::with_capacity(1 << k);
    for subset in 0..(1u8 << k) {
        let mut edited = codon;
        let mut edit_mask = 0u8;
        for (bit, &pos) in c_positions.iter().enumerate() {
            if subset & (1 << bit) != 0 {
                edited[pos] = b'T';
                edit_mask |= 1 << pos;
            }
        }
        out.push(CodonHypothesis {
            aa: translate_codon(edited),
            n_edits: subset.count_ones() as u8,
            edit_mask,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_standard_code() {
        assert_eq!(translate_codon(*b"ATG"), b'M');
        assert_eq!(translate_codon(*b"TAA"), STOP);
        assert_eq!(translate_codon(*b"CGA"), b'R');
    }

    #[test]
    fn acg_to_atg_start_is_one_of_the_hypotheses() {
        // ACG has no C-to-T-able C in the position that matters for this
        // classic case - the edit is at position 1 (A-C-G -> A-T-G).
        let hyps = codon_hypotheses(*b"ACG");
        assert!(hyps.iter().any(|h| h.aa == b'M' && h.n_edits == 1));
        assert!(hyps.iter().any(|h| h.aa == b'T' && h.n_edits == 0)); // unedited ACG = Thr
    }

    #[test]
    fn cga_can_resolve_to_edited_stop() {
        // CGA (Arg) -> TGA (stop) via a single C->T edit at position 0.
        let hyps = codon_hypotheses(*b"CGA");
        assert!(hyps.iter().any(|h| h.aa == STOP && h.n_edits == 1));
        assert!(hyps.iter().any(|h| h.aa == b'R' && h.n_edits == 0));
    }

    #[test]
    fn no_c_bases_has_a_single_hypothesis() {
        let hyps = codon_hypotheses(*b"ATG");
        assert_eq!(hyps.len(), 1);
        assert_eq!(hyps[0].n_edits, 0);
    }
}
