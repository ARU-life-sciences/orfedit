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

/// A frame change the nucleotide-level alignment (`align_nt`) made: `n_bases`
/// (1 or 2) sense bases starting at `sense_pos` belong to no codon, between
/// profile columns `after_column - 1` and `after_column`.
///
/// Skipping 1 base means the window has one base more than the reading
/// frame wants (a 1-bp insertion, or equivalently a 2-bp deletion); skipping
/// 2 means one base fewer (a 1-bp deletion, or a 2-bp insertion).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Frameshift {
    pub sense_pos: usize,
    pub n_bases: usize,
    pub after_column: usize,
}

/// One profile column's outcome in a nucleotide-level alignment: like
/// `ColumnCall`, but positioned by the codon's first sense base rather
/// than by a codon index in one frame.
#[derive(Debug, Clone)]
pub enum NtColumnCall {
    Aligned { sense_start: usize, aa: u8, n_edits: u8, edit_mask: u8 },
    Deleted,
}

pub struct NtAlignResult {
    pub score: f64,
    pub columns: Vec<NtColumnCall>,
    /// Sense-relative half-open [start, end) from the first aligned codon's
    /// first base to the last aligned codon's last base.
    pub sense_span: Option<(usize, usize)>,
    pub frameshifts: Vec<Frameshift>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum NtPtr {
    None,
    FromM,
    FromEx,
    FromEy,
    /// Ex entered from M by skipping 1 or 2 bases (a frameshift).
    Shift1,
    Shift2,
}

/// The same semi-global, edit-tolerant profile alignment as `align`, but over
/// the window's sense bases instead of one frame's codons, so the three
/// frames are one DP: a codon at sense position p is `sense_seq[p-3..p]`.
/// With `frameshift_penalty` infinite this finds the best single-frame
/// alignment (the best of `align` over the three frames - tested below).
/// With it finite, the alignment may also skip 1 or 2 bases between two
/// codons (from a match state, into the insertion state), paying that
/// penalty, which moves the rest of the gene into another frame - the way
/// a frameshifting indel would.
pub fn align_nt(
    profile: &Profile,
    sense_seq: &[u8],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: f64,
    frameshift_penalty: f64,
) -> NtAlignResult {
    let l = profile.len();
    let n = sense_seq.len();
    let w = n + 1;
    let idx = |i: usize, p: usize| i * w + p;
    let codon_at = |p: usize| [sense_seq[p - 3], sense_seq[p - 2], sense_seq[p - 1]];
    let shifts = frameshift_penalty.is_finite();

    // Scores roll over two profile rows; pointers are kept for traceback.
    let mut prev_m = vec![NEG_INF; w];
    let mut prev_ex = vec![0.0f64; w]; // free leading skip, any frame
    let mut prev_ey = vec![NEG_INF; w];
    for p in 0..w.min(3) {
        prev_m[p] = 0.0; // each frame's alignment may start at its own first base
    }
    let mut cur_m = vec![NEG_INF; w];
    let mut cur_ex = vec![NEG_INF; w];
    let mut cur_ey = vec![NEG_INF; w];
    let mut ptr_m = vec![NtPtr::None; (l + 1) * w];
    let mut ptr_ex = vec![NtPtr::None; (l + 1) * w];
    let mut ptr_ey = vec![NtPtr::None; (l + 1) * w];

    for i in 1..=l {
        for p in 0..w {
            // M[i][p]: column i-1 aligned to the codon ending at p.
            cur_m[p] = NEG_INF;
            if p >= 3 {
                let (pm, pex, pey) = (prev_m[p - 3], prev_ex[p - 3], prev_ey[p - 3]);
                let (best_prev, ptr) = if pm >= pex && pm >= pey {
                    (pm, NtPtr::FromM)
                } else if pex >= pey {
                    (pex, NtPtr::FromEx)
                } else {
                    (pey, NtPtr::FromEy)
                };
                if best_prev > NEG_INF {
                    let (hyp_score, ..) = best_hypothesis_score(profile, i - 1, codon_at(p), edit_penalty);
                    cur_m[p] = hyp_score + best_prev;
                    ptr_m[idx(i, p)] = ptr;
                }
            }
            // Ex[i][p]: an extra genome codon ending at p, or a frameshift.
            cur_ex[p] = NEG_INF;
            if p >= 3 {
                let open = cur_m[p - 3] - gap_open;
                let ext = cur_ex[p - 3] - gap_extend;
                if open >= ext {
                    cur_ex[p] = open;
                    ptr_ex[idx(i, p)] = NtPtr::FromM;
                } else {
                    cur_ex[p] = ext;
                    ptr_ex[idx(i, p)] = NtPtr::FromEx;
                }
            }
            if shifts && i < l {
                for (d, tag) in [(1usize, NtPtr::Shift1), (2, NtPtr::Shift2)] {
                    if p >= d {
                        let s = cur_m[p - d] - frameshift_penalty;
                        if s > cur_ex[p] {
                            cur_ex[p] = s;
                            ptr_ex[idx(i, p)] = tag;
                        }
                    }
                }
            }
            // Ey[i][p]: profile column i-1 with no genome codon.
            let open = prev_m[p] - gap_open;
            let ext = prev_ey[p] - gap_extend;
            if open >= ext {
                cur_ey[p] = open;
                ptr_ey[idx(i, p)] = NtPtr::FromM;
            } else {
                cur_ey[p] = ext;
                ptr_ey[idx(i, p)] = NtPtr::FromEy;
            }
        }
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_ex, &mut cur_ex);
        std::mem::swap(&mut prev_ey, &mut cur_ey);
    }

    if l == 0 {
        return NtAlignResult { score: 0.0, columns: vec![], sense_span: None, frameshifts: vec![] };
    }

    // Free trailing skip: prev_* now hold row L.
    let mut best_p = 0usize;
    let mut best_score = NEG_INF;
    let mut best_in_m = true;
    for p in 0..w {
        if prev_m[p] > best_score {
            best_score = prev_m[p];
            best_p = p;
            best_in_m = true;
        }
        if prev_ey[p] > best_score {
            best_score = prev_ey[p];
            best_p = p;
            best_in_m = false;
        }
    }

    let mut columns: Vec<NtColumnCall> = Vec::with_capacity(l);
    let mut frameshifts = Vec::new();
    let (mut i, mut p) = (l, best_p);
    let mut cur = if best_in_m { NtPtr::FromM } else { NtPtr::FromEy };
    let mut first: Option<usize> = None;
    let mut last: Option<usize> = None;
    while i > 0 {
        match cur {
            NtPtr::FromM => {
                let (_, aa, n_edits, edit_mask) = best_hypothesis_score(profile, i - 1, codon_at(p), edit_penalty);
                columns.push(NtColumnCall::Aligned { sense_start: p - 3, aa, n_edits, edit_mask });
                last.get_or_insert(p);
                first = Some(p - 3);
                cur = ptr_m[idx(i, p)];
                i -= 1;
                p -= 3;
            }
            NtPtr::FromEy => {
                columns.push(NtColumnCall::Deleted);
                cur = ptr_ey[idx(i, p)];
                i -= 1;
            }
            NtPtr::FromEx => match ptr_ex[idx(i, p)] {
                NtPtr::FromM => {
                    p -= 3;
                    cur = NtPtr::FromM;
                }
                NtPtr::FromEx => p -= 3,
                tag @ (NtPtr::Shift1 | NtPtr::Shift2) => {
                    let d = if tag == NtPtr::Shift1 { 1 } else { 2 };
                    frameshifts.push(Frameshift { sense_pos: p - d, n_bases: d, after_column: i });
                    p -= d;
                    cur = NtPtr::FromM;
                }
                _ => break,
            },
            _ => break, // the free leading skip (row 0)
        }
    }
    columns.reverse();
    frameshifts.reverse();
    let sense_span = match (first, last) {
        (Some(f), Some(la)) => Some((f, la)),
        _ => None,
    };
    NtAlignResult { score: best_score, columns, sense_span, frameshifts }
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

    /// Tiny deterministic generator, so the tests need no extra crate.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
    }

    fn random_dna(rng: &mut Lcg, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[(rng.next() % 4) as usize]).collect()
    }

    /// A random sense ORF (no internal stops) and the profile that wants it.
    fn random_gene(rng: &mut Lcg, n_codons: usize) -> (Vec<u8>, Profile) {
        let mut seq = Vec::new();
        let mut aas = Vec::new();
        while aas.len() < n_codons {
            let c = [b"ACGT"[(rng.next() % 4) as usize], b"ACGT"[(rng.next() % 4) as usize], b"ACGT"[(rng.next() % 4) as usize]];
            let aa = crate::codon::translate_codon(c);
            if aa != crate::codon::STOP {
                seq.extend_from_slice(&c);
                aas.push(aa);
            }
        }
        (seq, single_aa_profile("rand", &aas))
    }

    #[test]
    fn nt_alignment_without_shifts_is_the_best_single_frame() {
        let mut rng = Lcg(7);
        for _ in 0..40 {
            let n_codons = 20 + (rng.next() % 20) as usize;
            let (gene, mut profile) = random_gene(&mut rng, n_codons);
            // perturb the profile so it isn't a perfect match
            for col in profile.columns.iter_mut() {
                col[(rng.next() % 20) as usize] += 3.0;
            }
            let lead = (rng.next() % 40) as usize;
            let mut seq = random_dna(&mut rng, lead);
            seq.extend_from_slice(&gene);
            let trail = (rng.next() % 40) as usize;
            seq.extend(random_dna(&mut rng, trail));
            let nt = align_nt(&profile, &seq, 8.0, 1.0, 2.5, f64::INFINITY);
            let best = (0..3)
                .map(|f| {
                    let codons: Vec<[u8; 3]> = seq[f..].chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
                    let r = align(&profile, &codons, 8.0, 1.0, 2.5);
                    (r.score, r.codon_span.map(|(a, b)| (f + 3 * a, f + 3 * b)))
                })
                .fold((NEG_INF, None), |acc, x| if x.0 > acc.0 { x } else { acc });
            assert!((nt.score - best.0).abs() < 1e-9, "{} vs {}", nt.score, best.0);
            assert_eq!(nt.sense_span, best.1);
            assert!(nt.frameshifts.is_empty());
        }
    }

    #[test]
    fn finds_an_inserted_base() {
        let mut rng = Lcg(11);
        let (gene, profile) = random_gene(&mut rng, 60);
        let mut seq = random_dna(&mut rng, 30);
        let cut = seq.len() + 90; // after codon 30
        seq.extend_from_slice(&gene);
        seq.insert(cut, b'A');
        seq.extend(random_dna(&mut rng, 30));
        let r = align_nt(&profile, &seq, 8.0, 1.0, 2.5, 20.0);
        assert_eq!(r.frameshifts.len(), 1);
        let f = &r.frameshifts[0];
        assert_eq!(f.n_bases, 1);
        // the skipped base is at the insertion, give or take the codon it
        // lands in (equal-scoring placements within it)
        assert!(f.sense_pos.abs_diff(cut) <= 3, "{} vs {}", f.sense_pos, cut);
        assert_eq!(r.columns.iter().filter(|c| matches!(c, NtColumnCall::Deleted)).count(), 0);
        assert_eq!(r.sense_span, Some((30, 30 + gene.len() + 1)));
    }

    #[test]
    fn finds_a_deleted_base() {
        let mut rng = Lcg(13);
        let (gene, profile) = random_gene(&mut rng, 60);
        let mut seq = random_dna(&mut rng, 30);
        let cut = seq.len() + 120;
        seq.extend_from_slice(&gene);
        seq.remove(cut);
        seq.extend(random_dna(&mut rng, 30));
        let r = align_nt(&profile, &seq, 8.0, 1.0, 2.5, 20.0);
        assert_eq!(r.frameshifts.len(), 1);
        assert_eq!(r.frameshifts[0].n_bases, 2);
        assert!(r.frameshifts[0].sense_pos.abs_diff(cut) <= 3);
        // the codon broken by the deletion is the one profile column lost
        assert_eq!(r.columns.iter().filter(|c| matches!(c, NtColumnCall::Deleted)).count(), 1);
    }

    #[test]
    fn intact_genes_get_no_frameshifts() {
        let mut rng = Lcg(17);
        for _ in 0..40 {
            let (gene, profile) = random_gene(&mut rng, 50);
            let mut seq = random_dna(&mut rng, 50);
            seq.extend_from_slice(&gene);
            seq.extend(random_dna(&mut rng, 50));
            let r = align_nt(&profile, &seq, 8.0, 1.0, 2.5, 20.0);
            assert!(r.frameshifts.is_empty());
            assert_eq!(r.sense_span, Some((50, 50 + gene.len())));
        }
    }
}
