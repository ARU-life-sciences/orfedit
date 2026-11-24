//! Standard genetic code (NCBI translation table 1) plus C->U (C->T at the
//! DNA level) edit-hypothesis enumeration.
//!
//! Seed-plant mitochondrial/plastid RNA editing is (almost exclusively)
//! C->U at the transcript level, i.e. a genomic C can read as a T once
//! transcribed and edited. Ferns, lycophytes and hornworts also edit the
//! other way, U->C, and heavily: a Polypodiales mitochondrial gene has
//! about nine genomic in-frame stop codons that editing turns into Gln or
//! Arg. Those readings are a separate set (`codon_hypotheses_reverse`),
//! used only when the caller gives U->C edits a finite cost
//! ([`EditCosts`]). No other base change is ever hypothesized.

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
/// how many edits that reading assumes.
#[derive(Clone, Copy, Debug)]
pub struct CodonHypothesis {
    pub aa: u8,
    /// C->T (C->U) edits.
    pub n_edits: u8,
    /// T->C (U->C, "reverse") edits other than a stop read-through.
    pub n_rev: u8,
    /// 1 if this reading turns a genomic stop codon into a residue by a
    /// T->C edit of its first base (TAA, TAG -> Gln; TGA -> Arg).
    pub n_stop: u8,
    /// Bitmask (bits 0,1,2 = codon positions 0,1,2) of which positions were
    /// hypothesized edited, either way: a C there was read as T, a T as C.
    pub edit_mask: u8,
}

/// What an edit costs the alignment score. Reverse (U->C) edits come in
/// two kinds with a cost each: reading a genomic stop codon through
/// (`u_to_c_stop`) and any other (`u_to_c`). A stop inside a gene is
/// strong evidence by itself that an edit is due there - in fern
/// mitochondria nearly half of all U->C edits remove one - while an
/// ordinary T read as C is easily a chance improvement: at one flat cost
/// low enough to read stops through, seed-plant genes (which have no
/// reverse editing) picked up two false U->C edits each. Either cost
/// infinite (what a plain `f64` converts to, for both) rules that kind
/// out; both infinite is C->U editing only, as for seed plants and mosses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EditCosts {
    pub c_to_u: f64,
    pub u_to_c: f64,
    pub u_to_c_stop: f64,
    /// Costs of a C->U edit at a column of a profile that says where its
    /// reference proteins are edited (`Profile::edited`): at a column
    /// edited in them, and at one that is not. Both the same as `c_to_u`
    /// unless set apart; `c_to_u` is the cost at every column of a profile
    /// that does not say. Where a profile's proteins are RNA-edited ones,
    /// the columns their editing touches are known, and a C->U anywhere
    /// else is far more often a residue the species simply differs in: of
    /// the edits called in the plastid genes of 89 angiosperms at one flat
    /// cost, 42% were in their RNA; with a cost of 10 away from the known
    /// columns, 95%.
    pub c_to_u_known: f64,
    pub c_to_u_novel: f64,
}

impl EditCosts {
    /// These costs at a profile column: `site` is whether the references
    /// are edited there (`Profile::edit_site`), `None` for a profile that
    /// does not say, which leaves the plain cost.
    #[inline]
    pub fn at(self, site: Option<bool>) -> Self {
        match site {
            None => self,
            Some(true) => EditCosts { c_to_u: self.c_to_u_known, ..self },
            Some(false) => EditCosts { c_to_u: self.c_to_u_novel, ..self },
        }
    }
    /// Whether any reverse reading is allowed at all.
    pub fn reverse(&self) -> bool {
        self.u_to_c.is_finite() || self.u_to_c_stop.is_finite()
    }
    /// Cost of `n_edits` C->U edits, `n_rev` ordinary U->C edits and
    /// `n_stop` stop read-throughs; infinite if a kind that is ruled out
    /// is among them.
    #[inline]
    pub fn cost(&self, n_edits: u8, n_rev: u8, n_stop: u8) -> f64 {
        let mut cost = self.c_to_u * n_edits as f64;
        if n_rev > 0 {
            cost += self.u_to_c * n_rev as f64;
        }
        if n_stop > 0 {
            cost += self.u_to_c_stop * n_stop as f64;
        }
        cost
    }
}

impl From<f64> for EditCosts {
    fn from(c_to_u: f64) -> Self {
        EditCosts { c_to_u, u_to_c: f64::INFINITY, u_to_c_stop: f64::INFINITY, c_to_u_known: c_to_u, c_to_u_novel: c_to_u }
    }
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
            n_rev: 0,
            n_stop: 0,
            edit_mask,
        });
    }
    out
}

/// Every reading of `codon` with at least one T->C edit, each combined
/// with every subset of C->T edits at its C positions: the readings
/// `codon_hypotheses` leaves out. A stop codon's (TAA, TAG -> Gln; TGA ->
/// Arg) are among them. A reading never yields a stop here unless the
/// codon also needs a C->T edit to make one.
pub fn codon_hypotheses_reverse(codon: [u8; 3]) -> Vec<CodonHypothesis> {
    let ct_positions: Vec<usize> = (0..3).filter(|&p| codon[p] == b'C' || codon[p] == b'T').collect();
    let k = ct_positions.len();
    let mut out = Vec::new();
    for subset in 1..(1u8 << k) {
        let mut edited = codon;
        let (mut edit_mask, mut n_edits, mut n_rev) = (0u8, 0u8, 0u8);
        for (bit, &pos) in ct_positions.iter().enumerate() {
            if subset & (1 << bit) != 0 {
                if codon[pos] == b'C' {
                    edited[pos] = b'T';
                    n_edits += 1;
                } else {
                    edited[pos] = b'C';
                    n_rev += 1;
                }
                edit_mask |= 1 << pos;
            }
        }
        if n_rev > 0 {
            // a stop codon has one T, its first base, and no C: its only
            // reverse reading is the read-through
            let n_stop = (translate_codon(codon) == STOP) as u8;
            out.push(CodonHypothesis { aa: translate_codon(edited), n_edits, n_rev: n_rev - n_stop, n_stop, edit_mask });
        }
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
    fn reverse_readings_turn_stops_into_gln_and_arg() {
        for (codon, aa) in [(*b"TAA", b'Q'), (*b"TAG", b'Q'), (*b"TGA", b'R')] {
            let hyps = codon_hypotheses_reverse(codon);
            assert_eq!(hyps.len(), 1);
            assert_eq!((hyps[0].aa, hyps[0].n_edits, hyps[0].n_rev, hyps[0].n_stop, hyps[0].edit_mask), (aa, 0, 0, 1, 0b001));
            // and the forward set alone has only the stop
            assert!(codon_hypotheses(codon).iter().all(|h| h.aa == STOP && h.n_rev == 0));
        }
        // TCT: T->C at 0 or 2, alone, together, and each with C->T at 1
        let hyps = codon_hypotheses_reverse(*b"TCT");
        assert_eq!(hyps.len(), 6);
        assert!(hyps.iter().all(|h| h.n_rev >= 1 && h.n_stop == 0));
        assert!(hyps.iter().any(|h| h.aa == b'P' && h.n_rev == 1 && h.n_edits == 0)); // CCT
        assert!(hyps.iter().any(|h| h.aa == b'L' && h.n_rev == 1 && h.n_edits == 1)); // CTT
        // no T or C: nothing to add
        assert!(codon_hypotheses_reverse(*b"AAG").is_empty());
        assert!(codon_hypotheses_reverse(*b"CCA").is_empty());
    }

    #[test]
    fn edit_costs_from_a_plain_penalty_have_no_reverse_editing() {
        let costs = EditCosts::from(2.5);
        assert!(!costs.reverse());
        assert_eq!(costs.cost(2, 0, 0), 5.0);
        assert!(costs.cost(0, 0, 1).is_infinite() && costs.cost(0, 1, 0).is_infinite());
        let both = EditCosts { c_to_u: 2.5, u_to_c: 8.0, u_to_c_stop: 3.0, c_to_u_known: 2.5, c_to_u_novel: 2.5 };
        assert!(both.reverse());
        assert_eq!(both.cost(1, 2, 0), 18.5);
        assert_eq!(both.cost(0, 0, 1), 3.0);
        // stops only: an ordinary reverse edit stays ruled out
        let stops = EditCosts { c_to_u: 2.5, u_to_c: f64::INFINITY, u_to_c_stop: 3.0, c_to_u_known: 2.5, c_to_u_novel: 2.5 };
        // a profile that says where its references are edited: one cost for a C->U at such a
        // column, another elsewhere; nothing else changes, and nothing at all without the say
        let apart = EditCosts { c_to_u: 2.5, u_to_c: 8.0, u_to_c_stop: 3.0, c_to_u_known: 1.0, c_to_u_novel: 9.0 };
        assert_eq!(apart.at(None).cost(1, 1, 0), 10.5);
        assert_eq!(apart.at(Some(true)).cost(1, 1, 0), 9.0);
        assert_eq!(apart.at(Some(false)).cost(1, 1, 0), 17.0);
        assert!(stops.reverse());
        assert!(stops.cost(0, 1, 0).is_infinite());
        assert_eq!(stops.cost(0, 0, 1), 3.0);
    }

    #[test]
    fn no_c_bases_has_a_single_hypothesis() {
        let hyps = codon_hypotheses(*b"ATG");
        assert_eq!(hyps.len(), 1);
        assert_eq!(hyps[0].n_edits, 0);
    }
}
