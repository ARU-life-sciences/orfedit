//! Edit-tolerant profile-to-DNA semi-global alignment.
//!
//! Aligns a per-gene amino acid profile (global: the full profile must be
//! consumed) against one reading frame of a genomic window (semi-global:
//! leading/trailing genome codons outside the real gene are free - that's
//! what makes this a *boundary finder*, not just a scorer). Standard affine
//! gap DP, three matrices (M/Ex/Ey), except the match score at each cell is
//! the best-scoring C->T edit hypothesis for that codon (see `codon.rs`),
//! not just its literal genomic translation - that's what lets a genomic
//! sense codon resolve as an edited stop, or a genomic near-start resolve as
//! an edited ATG, whichever the reference profile supports better.

use crate::codon::codon_hypotheses;
use crate::profile::Profile;

const NEG_INF: f64 = f64::NEG_INFINITY;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Ptr {
    None,
    FromM,
    FromEx,
    FromEy,
}

/// One profile column's outcome after alignment.
#[derive(Debug, Clone)]
pub enum ColumnCall {
    /// Aligned to genome codon `codon_idx` (0-indexed within the codon
    /// sequence passed in), producing amino acid `aa` via `n_edits` C->T
    /// edits (`edit_mask`: bit p set = codon position p was hypothesized
    /// edited).
    Aligned { codon_idx: usize, aa: u8, n_edits: u8, edit_mask: u8 },
    /// Profile column has no corresponding genome codon (rare - a
    /// species-specific deletion relative to the reference profile).
    Deleted,
}

pub struct AlignResult {
    pub score: f64,
    /// One entry per profile column, in order.
    pub columns: Vec<ColumnCall>,
    /// 0-indexed [start, end) into the codon sequence spanned by aligned
    /// columns (i.e. the corrected gene body, in codons, within this frame's
    /// window). `None` if nothing aligned at all.
    pub codon_span: Option<(usize, usize)>,
}

fn best_hypothesis_score(profile: &Profile, col: usize, codon: [u8; 3], edit_penalty: f64) -> (f64, u8, u8, u8) {
    let mut best = (NEG_INF, b'X', 0u8, 0u8);
    for h in codon_hypotheses(codon) {
        let s = profile.score(col, h.aa) - edit_penalty * h.n_edits as f64;
        if s > best.0 {
            best = (s, h.aa, h.n_edits, h.edit_mask);
        }
    }
    best
}

/// `codons`: one reading frame's worth of complete codons from the genomic
/// window, in genomic order (index 0 = first/leftmost codon in the window).
pub fn align(
    profile: &Profile,
    codons: &[[u8; 3]],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: f64,
) -> AlignResult {
    let l = profile.len();
    let m = codons.len();
    let w = m + 1;
    let idx = |i: usize, j: usize| i * w + j;

    let mut mat_m = vec![NEG_INF; (l + 1) * w];
    let mut mat_ex = vec![NEG_INF; (l + 1) * w]; // extra genome codon (insertion rel. to profile)
    let mut mat_ey = vec![NEG_INF; (l + 1) * w]; // profile column with no genome codon (deletion)
    let mut ptr_m = vec![Ptr::None; (l + 1) * w];
    let mut ptr_ex = vec![Ptr::None; (l + 1) * w];
    let mut ptr_ey = vec![Ptr::None; (l + 1) * w];

    mat_m[idx(0, 0)] = 0.0;
    for j in 0..w {
        mat_ex[idx(0, j)] = 0.0; // free leading skip of genome codons before the gene starts
    }
    for i in 1..=l {
        // Ey column j=0: real deletion cost, profile side is NOT free-ended.
        let open = mat_m[idx(i - 1, 0)] - gap_open;
        let ext = mat_ey[idx(i - 1, 0)] - gap_extend;
        if open >= ext {
            mat_ey[idx(i, 0)] = open;
            ptr_ey[idx(i, 0)] = Ptr::FromM;
        } else {
            mat_ey[idx(i, 0)] = ext;
            ptr_ey[idx(i, 0)] = Ptr::FromEy;
        }
    }

    for i in 1..=l {
        for j in 1..=m {
            // M[i][j]
            let (hyp_score, aa, n_edits, edit_mask) =
                best_hypothesis_score(profile, i - 1, codons[j - 1], edit_penalty);
            let (pm, pex, pey) = (mat_m[idx(i - 1, j - 1)], mat_ex[idx(i - 1, j - 1)], mat_ey[idx(i - 1, j - 1)]);
            let (best_prev, ptr) = if pm >= pex && pm >= pey {
                (pm, Ptr::FromM)
            } else if pex >= pey {
                (pex, Ptr::FromEx)
            } else {
                (pey, Ptr::FromEy)
            };
            if best_prev > NEG_INF {
                mat_m[idx(i, j)] = hyp_score + best_prev;
                ptr_m[idx(i, j)] = ptr;
            }
            // stash the winning hypothesis inline via re-derivation at traceback time
            let _ = (aa, n_edits, edit_mask);

            // Ex[i][j]: consume a genome codon, profile index stays i.
            let open = mat_m[idx(i, j - 1)] - gap_open;
            let ext = mat_ex[idx(i, j - 1)] - gap_extend;
            if open >= ext {
                mat_ex[idx(i, j)] = open;
                ptr_ex[idx(i, j)] = Ptr::FromM;
            } else {
                mat_ex[idx(i, j)] = ext;
                ptr_ex[idx(i, j)] = Ptr::FromEx;
            }

            // Ey[i][j]: profile column i has no genome codon.
            let open = mat_m[idx(i - 1, j)] - gap_open;
            let ext = mat_ey[idx(i - 1, j)] - gap_extend;
            if open >= ext {
                mat_ey[idx(i, j)] = open;
                ptr_ey[idx(i, j)] = Ptr::FromM;
            } else {
                mat_ey[idx(i, j)] = ext;
                ptr_ey[idx(i, j)] = Ptr::FromEy;
            }
        }
    }

    if l == 0 {
        return AlignResult { score: 0.0, columns: vec![], codon_span: None };
    }

    // Free trailing skip: best end point is whichever j maximises M[L][j] or Ey[L][j].
    let mut best_j = 0usize;
    let mut best_score = NEG_INF;
    let mut best_in_m = true;
    for j in 0..=m {
        let sm = mat_m[idx(l, j)];
        if sm > best_score {
            best_score = sm;
            best_j = j;
            best_in_m = true;
        }
        let se = mat_ey[idx(l, j)];
        if se > best_score {
            best_score = se;
            best_j = j;
            best_in_m = false;
        }
    }

    // Traceback.
    let mut columns: Vec<ColumnCall> = Vec::with_capacity(l);
    let mut i = l;
    let mut j = best_j;
    let mut cur = if best_in_m { Ptr::FromM } else { Ptr::FromEy }; // "current matrix", reusing Ptr as a 2-value tag
    let mut first_codon_idx: Option<usize> = None;
    let mut last_codon_idx: Option<usize> = None;

    while i > 0 {
        match cur {
            Ptr::FromM => {
                let (_, aa, n_edits, edit_mask) =
                    best_hypothesis_score(profile, i - 1, codons[j - 1], edit_penalty);
                columns.push(ColumnCall::Aligned { codon_idx: j - 1, aa, n_edits, edit_mask });
                last_codon_idx.get_or_insert(j - 1);
                first_codon_idx = Some(j - 1);
                let p = ptr_m[idx(i, j)];
                i -= 1;
                j -= 1;
                cur = p;
            }
            Ptr::FromEy => {
                columns.push(ColumnCall::Deleted);
                let p = ptr_ey[idx(i, j)];
                i -= 1;
                cur = p;
            }
            Ptr::FromEx => {
                let p = ptr_ex[idx(i, j)];
                j -= 1;
                cur = p;
            }
            Ptr::None => break, // reached the free-skip zone (i==0 boundary) or malformed path
        }
    }
    columns.reverse();

    let codon_span = match (first_codon_idx, last_codon_idx) {
        (Some(f), Some(la)) => Some((f, la + 1)),
        _ => None,
    };

    AlignResult { score: best_score, columns, codon_span }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn codons_from_str(s: &str) -> Vec<[u8; 3]> {
        let b = s.as_bytes();
        b.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
    }

    fn single_aa_profile(gene: &str, aas: &[u8]) -> Profile {
        use crate::codon::AA_ORDER;
        let columns = aas
            .iter()
            .map(|&target| {
                let mut col = [-5.0f64; 20];
                let idx = AA_ORDER.iter().position(|&x| x == target).unwrap();
                col[idx] = 5.0;
                col
            })
            .collect();
        Profile { gene: gene.to_string(), columns }
    }

    #[test]
    fn rescues_acg_to_atg_start() {
        // Profile expects M-K-V. Genome has flanking junk, then ACG-AAA-GTT
        // (ACG = unedited Thr, needs a 1-edit read as Met to match).
        let profile = single_aa_profile("test", b"MKV");
        let genome = "TTTTTTACGAAAGTTCCCCCC"; // junk(2 codons) + ACG AAA GTT + junk(2 codons)
        let codons = codons_from_str(genome);
        let r = align(&profile, &codons, 4.0, 1.0, 2.0);
        assert_eq!(r.columns.len(), 3);
        match &r.columns[0] {
            ColumnCall::Aligned { aa, n_edits, codon_idx, .. } => {
                assert_eq!(*aa, b'M');
                assert_eq!(*n_edits, 1);
                assert_eq!(*codon_idx, 2); // 3rd codon in the window (0-indexed 2)
            }
            other => panic!("expected aligned M, got {:?}", other),
        }
    }

    #[test]
    fn resolves_edited_internal_stop_as_sense_codon() {
        // Profile expects M-R-V. Genome mid-codon is CGA (Arg, unedited) -
        // must NOT be reinterpreted as a stop, since profile wants Arg here
        // and Arg is CGA's unedited reading (0 edits beats any stop reading).
        let profile = single_aa_profile("test", b"MRV");
        let genome = "ATGCGAGTT";
        let codons = codons_from_str(genome);
        let r = align(&profile, &codons, 4.0, 1.0, 2.0);
        assert_eq!(r.columns.len(), 3);
        match &r.columns[1] {
            ColumnCall::Aligned { aa, n_edits, .. } => {
                assert_eq!(*aa, b'R');
                assert_eq!(*n_edits, 0);
            }
            other => panic!("expected aligned R, got {:?}", other),
        }
    }

    #[test]
    fn unrescuable_interior_stop_is_deleted_not_aligned() {
        // Profile expects M-K-V. Genome has ATG(M, exact) then TAA - a
        // "pure" stop with zero Cs, so no edit hypothesis can rescue it -
        // then GTT(V, exact). The DP must explain the TAA codon as a
        // skipped/deleted position, never as an aligned "K" - a real
        // regression: an earlier version scored unrecognised amino acids
        // (stop, X) at a finite -10.0, which could let a stop win an
        // interior column over a deletion and leak a literal '*' into
        // output protein sequences (see profile.rs).
        let profile = single_aa_profile("test", b"MKV");
        let codons = codons_from_str("ATGTAAGTT");
        let r = align(&profile, &codons, 8.0, 1.0, 2.5);
        assert_eq!(r.columns.len(), 3);
        match &r.columns[1] {
            ColumnCall::Deleted => {}
            other => panic!("expected column 1 (K) to be Deleted, got {:?}", other),
        }
        // And no Aligned column should ever carry a stop/unknown character.
        for col in &r.columns {
            if let ColumnCall::Aligned { aa, .. } = col {
                assert_ne!(*aa, crate::codon::STOP);
                assert_ne!(*aa, crate::codon::UNKNOWN);
            }
        }
    }

    #[test]
    fn no_change_case_is_stable() {
        let profile = single_aa_profile("test", b"MKV");
        let codons = codons_from_str("ATGAAAGTT");
        let r = align(&profile, &codons, 4.0, 1.0, 2.0);
        assert_eq!(r.codon_span, Some((0, 3)));
        for col in &r.columns {
            match col {
                ColumnCall::Aligned { n_edits, .. } => assert_eq!(*n_edits, 0),
                ColumnCall::Deleted => panic!("unexpected deletion in exact-match case"),
            }
        }
    }
}
