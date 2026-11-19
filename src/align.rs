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

use crate::codon::{codon_hypotheses, codon_hypotheses_reverse, EditCosts};
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

/// The best reading of `codon` against profile column `col`: (score, amino
/// acid, number of edits either way, mask of edited positions). C->U
/// readings first, so that without reverse editing - and on a tie with it -
/// the result is what it always was.
fn best_hypothesis_score(profile: &Profile, col: usize, codon: [u8; 3], edit_penalty: EditCosts) -> (f64, u8, u8, u8) {
    let mut best = (NEG_INF, b'X', 0u8, 0u8);
    for h in codon_hypotheses(codon) {
        let s = profile.score(col, h.aa) - edit_penalty.cost(h.n_edits, 0, 0);
        if s > best.0 {
            best = (s, h.aa, h.n_edits, h.edit_mask);
        }
    }
    if edit_penalty.reverse() {
        for h in codon_hypotheses_reverse(codon) {
            let s = profile.score(col, h.aa) - edit_penalty.cost(h.n_edits, h.n_rev, h.n_stop);
            if s > best.0 {
                best = (s, h.aa, h.n_edits + h.n_rev + h.n_stop, h.edit_mask);
            }
        }
    }
    best
}

/// Base class for the codon-score table: A, C, G, T, anything else.
/// `translate_codon` and `codon_hypotheses` treat every non-ACGT byte the
/// same way (it only ever matches a wildcard), so five classes are exact.
#[inline]
fn base_class(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}
const N_CODON_CLASSES: usize = 125;

/// For each of the 125 codon classes, its edit hypotheses as
/// (index into AA_ORDER or None, C->U edits, U->C edits, stop
/// read-throughs): the C->U readings in `codon_hypotheses` order, then the
/// reverse ones.
type ClassReading = (Option<usize>, u8, u8, u8);
type ClassHypotheses = (Vec<ClassReading>, Vec<ClassReading>);
fn codon_class_hypotheses() -> &'static Vec<ClassHypotheses> {
    static TABLE: std::sync::OnceLock<Vec<ClassHypotheses>> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        const REP: [u8; 5] = *b"ACGTN";
        (0..N_CODON_CLASSES)
            .map(|c| {
                let codon = [REP[c / 25], REP[(c / 5) % 5], REP[c % 5]];
                let pack = |hyps: Vec<crate::codon::CodonHypothesis>| {
                    hyps.into_iter().map(|h| (crate::codon::aa_index(h.aa), h.n_edits, h.n_rev, h.n_stop)).collect::<Vec<_>>()
                };
                (pack(codon_hypotheses(codon)), pack(codon_hypotheses_reverse(codon)))
            })
            .collect()
    })
}

/// `best_hypothesis_score(profile, col, codon, edit_penalty).0` for every
/// profile column and codon class: `table[col * 125 + class]`. Same
/// arithmetic in the same order as `best_hypothesis_score`, so the values
/// are bit-identical; the DP then does one lookup per cell instead of
/// enumerating (and heap-allocating) the hypotheses each time.
fn hypothesis_score_table(profile: &Profile, edit_penalty: EditCosts) -> Vec<f64> {
    let hyps = codon_class_hypotheses();
    let reverse = edit_penalty.reverse();
    let mut table = Vec::with_capacity(profile.len() * N_CODON_CLASSES);
    for col in &profile.columns {
        for (forward, backward) in hyps {
            let mut best = NEG_INF;
            let readings = forward.iter().chain(backward.iter().filter(|_| reverse));
            for &(aa_idx, n_edits, n_rev, n_stop) in readings {
                let base = match aa_idx {
                    Some(a) => col[a],
                    None => NEG_INF,
                };
                let s = base - edit_penalty.cost(n_edits, n_rev, n_stop);
                if s > best {
                    best = s;
                }
            }
            table.push(best);
        }
    }
    table
}

/// Codon class of the codon ending at each sense position p (p >= 3).
fn codon_classes(sense_seq: &[u8]) -> Vec<u8> {
    let mut out = vec![0u8; sense_seq.len() + 1];
    for p in 3..=sense_seq.len() {
        out[p] = (base_class(sense_seq[p - 3]) * 25 + base_class(sense_seq[p - 2]) * 5 + base_class(sense_seq[p - 1])) as u8;
    }
    out
}

/// `codons`: one reading frame's worth of complete codons from the genomic
/// window, in genomic order (index 0 = first/leftmost codon in the window).
pub fn align(
    profile: &Profile,
    codons: &[[u8; 3]],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: impl Into<EditCosts>,
) -> AlignResult {
    let edit_penalty: EditCosts = edit_penalty.into();
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
    /// Sense-relative half-open [start, end) over which the frame change
    /// could equally well sit: moving it to any column boundary within
    /// the ungapped stretch either side (FRAMESHIFT_WINDOW_COLUMNS) that
    /// scores within FRAMESHIFT_WINDOW_TOLERANCE of the best placement.
    /// Where both frames read about as well as each other (a divergent
    /// stretch), that's wider than the bases skipped.
    pub window: (usize, usize),
}

/// How many aligned columns either side of a frame change are tried as
/// its other possible placements, and how close (score) to the best one
/// such a placement must come to widen the window.
pub const FRAMESHIFT_WINDOW_COLUMNS: usize = 30;
pub const FRAMESHIFT_WINDOW_TOLERANCE: f64 = 2.0;

/// The window (see `Frameshift::window`) of the frame change between
/// columns `c - 1` and `c`: within the ungapped run of consecutive codons
/// on each side, score every breakpoint k (columns before k read in the
/// upstream frame, from k on in the downstream one).
fn frameshift_window(
    profile: &Profile,
    sense_seq: &[u8],
    columns: &[NtColumnCall],
    c: usize,
    edit_penalty: impl Into<EditCosts>,
) -> Option<(usize, usize)> {
    let edit_penalty: EditCosts = edit_penalty.into();
    let start_of = |j: usize| match columns.get(j) {
        Some(NtColumnCall::Aligned { sense_start, .. }) => Some(*sense_start),
        _ => None,
    };
    let (a0, b0) = (start_of(c - 1)?, start_of(c)?); // the codons either side
    // the run of consecutive codons in each frame
    let mut lo = c - 1;
    while lo > 0 && c - lo < FRAMESHIFT_WINDOW_COLUMNS && start_of(lo - 1).is_some_and(|s| s + 3 == start_of(lo).unwrap()) {
        lo -= 1;
    }
    let mut hi = c;
    while hi + 1 < columns.len()
        && hi + 1 - c < FRAMESHIFT_WINDOW_COLUMNS
        && start_of(hi + 1).is_some_and(|s| s == start_of(hi).unwrap() + 3)
    {
        hi += 1;
    }
    // codon start of column j read in frame A (upstream) or B (downstream)
    let in_a = |j: usize| (a0 as i64 + 3 * (j as i64 - (c as i64 - 1))) as usize;
    let in_b = |j: usize| b0 as i64 + 3 * (j as i64 - c as i64);
    let score = |j: usize, s: i64| -> f64 {
        if s < 0 || s as usize + 3 > sense_seq.len() {
            return NEG_INF;
        }
        let s = s as usize;
        best_hypothesis_score(profile, j, [sense_seq[s], sense_seq[s + 1], sense_seq[s + 2]], edit_penalty).0
    };
    // breakpoints k in lo+1..=hi: columns lo..k in A, k..=hi in B
    let mut scores = Vec::new();
    for k in lo + 1..=hi {
        let a: f64 = (lo..k).map(|j| score(j, in_a(j) as i64)).sum();
        let b: f64 = (k..=hi).map(|j| score(j, in_b(j))).sum();
        scores.push((k, a + b));
    }
    let best = scores.iter().map(|x| x.1).fold(NEG_INF, f64::max);
    if !best.is_finite() {
        return None;
    }
    let near: Vec<usize> = scores.iter().filter(|x| x.1 >= best - FRAMESHIFT_WINDOW_TOLERANCE).map(|x| x.0).collect();
    // bases between the last upstream codon and the first downstream one
    let start = near.iter().map(|&k| in_a(k - 1) + 3).min()?;
    let end = near.iter().map(|&k| in_b(k).max(0) as usize).max()?;
    Some((start.min(end), end.max(start)))
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

#[cfg(test)]
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
    edit_penalty: impl Into<EditCosts>,
    frameshift_penalty: f64,
) -> NtAlignResult {
    let edit_penalty: EditCosts = edit_penalty.into();
    align_nt_from(profile, sense_seq, gap_open, gap_extend, edit_penalty, frameshift_penalty, true)
}

/// `align_nt` with the alignment starting in frame 0 of `sense_seq` (a
/// sequence already in frame, such as joined exons): with
/// `frameshift_penalty` infinite this is `align` over frame 0's codons.
pub fn align_nt_frame0(
    profile: &Profile,
    sense_seq: &[u8],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: impl Into<EditCosts>,
    frameshift_penalty: f64,
) -> NtAlignResult {
    let edit_penalty: EditCosts = edit_penalty.into();
    align_nt_from(profile, sense_seq, gap_open, gap_extend, edit_penalty, frameshift_penalty, false)
}

const PTR_M: u8 = 1;
const PTR_EX: u8 = 2;
const PTR_EY: u8 = 3;
const PTR_SHIFT1: u8 = 4;
const PTR_SHIFT2: u8 = 5;

/// Alignments with at least this many DP cells (profile columns x sense
/// bases) are done in two passes - see `align_nt_impl`.
const TWO_PASS_MIN_CELLS: usize = 8_000_000;

fn align_nt_from(
    profile: &Profile,
    sense_seq: &[u8],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: impl Into<EditCosts>,
    frameshift_penalty: f64,
    any_frame: bool,
) -> NtAlignResult {
    let edit_penalty: EditCosts = edit_penalty.into();
    align_nt_impl(profile, sense_seq, gap_open, gap_extend, edit_penalty, frameshift_penalty, any_frame, TWO_PASS_MIN_CELLS)
}

/// Where the best alignment lies, found without keeping any traceback.
struct Located {
    score: f64,
    /// Sense position the alignment ends at, and whether in M (else Ey).
    end_p: usize,
    in_m: bool,
    /// Sense position (in row 0) the alignment's traceback would end at.
    origin: usize,
}

/// One DP column (a sense position) of `locate`: the three states for
/// every profile row, and the origin each cell's path started from.
#[derive(Default)]
struct LocateColumn {
    m: Vec<f64>,
    ex: Vec<f64>,
    ey: Vec<f64>,
    om: Vec<u32>,
    oex: Vec<u32>,
    oey: Vec<u32>,
}

/// The DP of `align_nt_impl` over the whole of `sense_seq`, scores only:
/// each cell also carries the row-0 position its path started from (chosen
/// by exactly the comparisons that set the traceback pointers). Cells are
/// filled a sense position at a time (down the profile), which needs only
/// the previous three positions: four short columns that stay in cache,
/// however long the sequence, instead of rows as long as the contig.
/// Every cell is computed from the same operands by the same operations as
/// in the row-by-row fill, so scores are bit-identical.
#[allow(clippy::too_many_arguments)]
fn locate(
    l: usize,
    n: usize,
    hyp_table: &[f64],
    classes: &[u8],
    gap_open: f64,
    gap_extend: f64,
    frameshift_penalty: f64,
    any_frame: bool,
) -> Located {
    let w = n + 1;
    assert!(w < u32::MAX as usize, "sequence too long");
    let shifts = frameshift_penalty.is_finite();
    let m_start = if any_frame { 3 } else { 1 };
    // hypothesis scores by codon class, then profile column: one contiguous
    // run per sense position
    let mut hyp_by_class = vec![0.0f64; N_CODON_CLASSES * l];
    for col in 0..l {
        for class in 0..N_CODON_CLASSES {
            hyp_by_class[class * l + col] = hyp_table[col * N_CODON_CLASSES + class];
        }
    }
    let mut ring: Vec<LocateColumn> = (0..4)
        .map(|_| LocateColumn {
            m: vec![NEG_INF; l + 1],
            ex: vec![NEG_INF; l + 1],
            ey: vec![NEG_INF; l + 1],
            om: vec![0; l + 1],
            oex: vec![0; l + 1],
            oey: vec![0; l + 1],
        })
        .collect();
    let mut loc = Located { score: NEG_INF, end_p: 0, in_m: true, origin: 0 };
    // In one frame without frame changes, only every third position can be
    // reached at all (the rest stay -inf throughout): skip the others.
    let step = if !any_frame && !shifts { 3 } else { 1 };
    for p in (0..w).step_by(step) {
        let mut c = std::mem::take(&mut ring[p % 4]);
        // row 0: the free leading skip, and each frame's own first base
        c.m[0] = if p < m_start { 0.0 } else { NEG_INF };
        c.om[0] = p as u32;
        c.ex[0] = if any_frame || p % 3 == 0 { 0.0 } else { NEG_INF };
        c.oex[0] = p as u32;
        c.ey[0] = NEG_INF;
        if p >= 3 {
            let b = &ring[(p + 1) % 4]; // position p - 3
            let class = classes[p] as usize;
            let hyp = &hyp_by_class[class * l..][..l];
            // (Plain branches: measured faster here than branch-free
            // selects. An unreachable predecessor (-inf) gives -inf without
            // a test, since no hypothesis score is +inf or NaN.)
            {
                let (bm, bex, bey) = (&b.m[..l], &b.ex[..l], &b.ey[..l]);
                let (bom, boex, boey) = (&b.om[..l], &b.oex[..l], &b.oey[..l]);
                let (cm, com) = (&mut c.m[1..][..l], &mut c.om[1..][..l]);
                for i in 0..l {
                    let (pm, pex, pey) = (bm[i], bex[i], bey[i]);
                    let (best_prev, origin) = if pm >= pex && pm >= pey {
                        (pm, bom[i])
                    } else if pex >= pey {
                        (pex, boex[i])
                    } else {
                        (pey, boey[i])
                    };
                    cm[i] = hyp[i] + best_prev;
                    com[i] = origin;
                }
            }
            {
                let (bm, bex) = (&b.m[1..][..l], &b.ex[1..][..l]);
                let (bom, boex) = (&b.om[1..][..l], &b.oex[1..][..l]);
                let (cex, coex) = (&mut c.ex[1..][..l], &mut c.oex[1..][..l]);
                for i in 0..l {
                    let open = bm[i] - gap_open;
                    let ext = bex[i] - gap_extend;
                    if open >= ext {
                        cex[i] = open;
                        coex[i] = bom[i];
                    } else {
                        cex[i] = ext;
                        coex[i] = boex[i];
                    }
                }
            }
        } else {
            c.m[1..].fill(NEG_INF);
            c.ex[1..].fill(NEG_INF);
        }
        if shifts {
            for d in [1usize, 2] {
                if p >= d {
                    let b = &ring[(p + 4 - d) % 4]; // position p - d
                    for i in 1..l {
                        let s = b.m[i] - frameshift_penalty;
                        if s > c.ex[i] {
                            c.ex[i] = s;
                            c.oex[i] = b.om[i];
                        }
                    }
                }
            }
        }
        {
            let (mut ey, mut ey_origin) = (c.ey[0], c.oey[0]);
            let (cm, com) = (&c.m[..l], &c.om[..l]);
            let (cey, coey) = (&mut c.ey[1..][..l], &mut c.oey[1..][..l]);
            for i in 0..l {
                let open = cm[i] - gap_open;
                let ext = ey - gap_extend;
                if open >= ext {
                    ey = open;
                    ey_origin = com[i];
                } else {
                    ey = ext;
                }
                cey[i] = ey;
                coey[i] = ey_origin;
            }
        }
        if c.m[l] > loc.score {
            loc = Located { score: c.m[l], end_p: p, in_m: true, origin: c.om[l] as usize };
        }
        if c.ey[l] > loc.score {
            loc = Located { score: c.ey[l], end_p: p, in_m: false, origin: c.oey[l] as usize };
        }
        ring[p % 4] = c;
    }
    loc
}

/// The alignment DP proper. Small problems fill the whole matrix, keeping a
/// traceback pointer per cell. Large ones (`two_pass_min_cells`; whole-contig
/// scans) first run `locate` to learn where the best alignment starts and
/// ends, then fill only the columns between: the best path lies inside that
/// window, every other cell there scores no more than it does in the full
/// matrix, and all the DP's comparisons keep the on-path choice on a tie, so
/// the traceback - and hence the result - is the same as the full matrix's
/// (tested against the one-pass reference below).
#[allow(clippy::too_many_arguments)]
fn align_nt_impl(
    profile: &Profile,
    sense_seq: &[u8],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: impl Into<EditCosts>,
    frameshift_penalty: f64,
    any_frame: bool,
    two_pass_min_cells: usize,
) -> NtAlignResult {
    let edit_penalty: EditCosts = edit_penalty.into();
    let l = profile.len();
    let n = sense_seq.len();
    if l == 0 {
        return NtAlignResult { score: 0.0, columns: vec![], sense_span: None, frameshifts: vec![] };
    }
    let hyp_table = hypothesis_score_table(profile, edit_penalty);
    let classes = codon_classes(sense_seq);
    let located = if (l + 1).saturating_mul(n + 1) >= two_pass_min_cells {
        Some(locate(l, n, &hyp_table, &classes, gap_open, gap_extend, frameshift_penalty, any_frame))
    } else {
        None
    };
    // Sense positions lo..=hi are the DP columns filled (local q = p - lo).
    let (lo, hi) = located.as_ref().map_or((0, n), |loc| (loc.origin, loc.end_p));
    let w = hi - lo + 1;
    let idx = |i: usize, q: usize| i * w + q;
    let codon_at = |p: usize| [sense_seq[p - 3], sense_seq[p - 2], sense_seq[p - 1]];
    let shifts = frameshift_penalty.is_finite();
    let m_start = if any_frame { 3 } else { 1 };

    // Scores roll over two profile rows; pointers are kept for traceback.
    let mut prev_m = vec![NEG_INF; w];
    let mut prev_ex = vec![NEG_INF; w];
    let mut prev_ey = vec![NEG_INF; w];
    for q in 0..w {
        let p = lo + q;
        if any_frame || p % 3 == 0 {
            prev_ex[q] = 0.0; // free leading skip
        }
        if p < m_start {
            prev_m[q] = 0.0; // each frame's alignment may start at its own first base
        }
    }
    let mut cur_m = vec![NEG_INF; w];
    let mut cur_ex = vec![NEG_INF; w];
    let mut cur_ey = vec![NEG_INF; w];
    // Traceback pointers, one byte per cell: bits 0-1 the M pointer, bits
    // 2-4 the Ex pointer, bits 5-6 the Ey pointer (PTR_* codes below).
    let mut ptrs = vec![0u8; (l + 1) * w];
    let classes = &classes[lo..=hi];

    for i in 1..=l {
        let hyp_row = &hyp_table[(i - 1) * N_CODON_CLASSES..i * N_CODON_CLASSES];
        let ptr_row = &mut ptrs[i * w..(i + 1) * w];
        for p in 0..w {
            let mut cell = 0u8;
            // M[i][p]: column i-1 aligned to the codon ending at p.
            cur_m[p] = NEG_INF;
            if p >= 3 {
                let (pm, pex, pey) = (prev_m[p - 3], prev_ex[p - 3], prev_ey[p - 3]);
                let (best_prev, ptr) = if pm >= pex && pm >= pey {
                    (pm, PTR_M)
                } else if pex >= pey {
                    (pex, PTR_EX)
                } else {
                    (pey, PTR_EY)
                };
                if best_prev > NEG_INF {
                    let hyp_score = hyp_row[classes[p] as usize];
                    cur_m[p] = hyp_score + best_prev;
                    cell |= ptr;
                }
            }
            // Ex[i][p]: an extra genome codon ending at p, or a frameshift.
            cur_ex[p] = NEG_INF;
            let mut ex_ptr = 0u8;
            if p >= 3 {
                let open = cur_m[p - 3] - gap_open;
                let ext = cur_ex[p - 3] - gap_extend;
                if open >= ext {
                    cur_ex[p] = open;
                    ex_ptr = PTR_M;
                } else {
                    cur_ex[p] = ext;
                    ex_ptr = PTR_EX;
                }
            }
            if shifts && i < l {
                for (d, tag) in [(1usize, PTR_SHIFT1), (2, PTR_SHIFT2)] {
                    if p >= d {
                        let s = cur_m[p - d] - frameshift_penalty;
                        if s > cur_ex[p] {
                            cur_ex[p] = s;
                            ex_ptr = tag;
                        }
                    }
                }
            }
            cell |= ex_ptr << 2;
            // Ey[i][p]: profile column i-1 with no genome codon.
            let open = prev_m[p] - gap_open;
            let ext = prev_ey[p] - gap_extend;
            if open >= ext {
                cur_ey[p] = open;
                cell |= PTR_M << 5;
            } else {
                cur_ey[p] = ext;
                cell |= PTR_EY << 5;
            }
            ptr_row[p] = cell;
        }
        std::mem::swap(&mut prev_m, &mut cur_m);
        std::mem::swap(&mut prev_ex, &mut cur_ex);
        std::mem::swap(&mut prev_ey, &mut cur_ey);
    }

    // Free trailing skip: prev_* now hold row L.
    let mut best_q = 0usize;
    let mut best_score = NEG_INF;
    let mut best_in_m = true;
    for q in 0..w {
        if prev_m[q] > best_score {
            best_score = prev_m[q];
            best_q = q;
            best_in_m = true;
        }
        if prev_ey[q] > best_score {
            best_score = prev_ey[q];
            best_q = q;
            best_in_m = false;
        }
    }
    if let Some(loc) = &located {
        // The window must reproduce the located optimum exactly.
        assert!(
            best_score.to_bits() == loc.score.to_bits() && lo + best_q == loc.end_p && best_in_m == loc.in_m,
            "two-pass alignment disagrees with its locating pass"
        );
    }

    let mut columns: Vec<NtColumnCall> = Vec::with_capacity(l);
    let mut frameshifts = Vec::new();
    let (mut i, mut q) = (l, best_q);
    let ptr_m = |i: usize, q: usize| ptrs[idx(i, q)] & 3;
    let ptr_ex = |i: usize, q: usize| (ptrs[idx(i, q)] >> 2) & 7;
    let ptr_ey = |i: usize, q: usize| (ptrs[idx(i, q)] >> 5) & 3;
    let mut cur = if best_in_m { PTR_M } else { PTR_EY };
    let mut first: Option<usize> = None;
    let mut last: Option<usize> = None;
    while i > 0 {
        let p = lo + q; // sense position
        match cur {
            PTR_M => {
                let (_, aa, n_edits, edit_mask) = best_hypothesis_score(profile, i - 1, codon_at(p), edit_penalty);
                columns.push(NtColumnCall::Aligned { sense_start: p - 3, aa, n_edits, edit_mask });
                last.get_or_insert(p);
                first = Some(p - 3);
                cur = ptr_m(i, q);
                i -= 1;
                q -= 3;
            }
            PTR_EY => {
                columns.push(NtColumnCall::Deleted);
                cur = ptr_ey(i, q);
                i -= 1;
            }
            PTR_EX => match ptr_ex(i, q) {
                PTR_M => {
                    q -= 3;
                    cur = PTR_M;
                }
                PTR_EX => q -= 3,
                tag @ (PTR_SHIFT1 | PTR_SHIFT2) => {
                    let d = if tag == PTR_SHIFT1 { 1 } else { 2 };
                    frameshifts.push(Frameshift { sense_pos: p - d, n_bases: d, after_column: i, window: (p - d, p) });
                    q -= d;
                    cur = PTR_M;
                }
                _ => break,
            },
            _ => break, // the free leading skip (row 0)
        }
    }
    columns.reverse();
    frameshifts.reverse();
    for f in frameshifts.iter_mut() {
        if let Some((s, e)) = frameshift_window(profile, sense_seq, &columns, f.after_column, edit_penalty) {
            f.window = (s.min(f.sense_pos), e.max(f.sense_pos + f.n_bases));
        }
    }
    let sense_span = match (first, last) {
        (Some(f), Some(la)) => Some((f, la)),
        _ => None,
    };
    NtAlignResult { score: best_score, columns, sense_span, frameshifts }
}

/// The score `align_nt_frame0` would return, without building the
/// alignment (no traceback matrix): for callers that only compare scores.
pub fn align_nt_frame0_score(
    profile: &Profile,
    sense_seq: &[u8],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: impl Into<EditCosts>,
    frameshift_penalty: f64,
) -> f64 {
    let edit_penalty: EditCosts = edit_penalty.into();
    if profile.is_empty() {
        return 0.0;
    }
    let hyp_table = hypothesis_score_table(profile, edit_penalty);
    let classes = codon_classes(sense_seq);
    locate(profile.len(), sense_seq.len(), &hyp_table, &classes, gap_open, gap_extend, frameshift_penalty, false).score
}

#[cfg(test)]
fn align_nt_from_reference(
    profile: &Profile,
    sense_seq: &[u8],
    gap_open: f64,
    gap_extend: f64,
    edit_penalty: impl Into<EditCosts>,
    frameshift_penalty: f64,
    any_frame: bool,
) -> NtAlignResult {
    let edit_penalty: EditCosts = edit_penalty.into();
    let l = profile.len();
    let n = sense_seq.len();
    let w = n + 1;
    let idx = |i: usize, p: usize| i * w + p;
    let codon_at = |p: usize| [sense_seq[p - 3], sense_seq[p - 2], sense_seq[p - 1]];
    let shifts = frameshift_penalty.is_finite();

    // Scores roll over two profile rows; pointers are kept for traceback.
    let mut prev_m = vec![NEG_INF; w];
    let mut prev_ex = vec![NEG_INF; w];
    let mut prev_ey = vec![NEG_INF; w];
    for p in 0..w {
        if any_frame || p % 3 == 0 {
            prev_ex[p] = 0.0; // free leading skip
        }
    }
    for p in 0..w.min(if any_frame { 3 } else { 1 }) {
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
                    frameshifts.push(Frameshift { sense_pos: p - d, n_bases: d, after_column: i, window: (p - d, p) });
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
    for f in frameshifts.iter_mut() {
        if let Some((s, e)) = frameshift_window(profile, sense_seq, &columns, f.after_column, edit_penalty) {
            f.window = (s.min(f.sense_pos), e.max(f.sense_pos + f.n_bases));
        }
    }
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

    #[test]
    fn frame0_alignment_without_shifts_is_align_over_frame0() {
        let mut rng = Lcg(23);
        for _ in 0..30 {
            let n_codons = 20 + (rng.next() % 20) as usize;
            let (gene, mut profile) = random_gene(&mut rng, n_codons);
            for col in profile.columns.iter_mut() {
                col[(rng.next() % 20) as usize] += 3.0;
            }
            let lead = (rng.next() % 40) as usize;
            let mut seq = random_dna(&mut rng, lead);
            seq.extend_from_slice(&gene);
            let codons: Vec<[u8; 3]> = seq.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
            let r = align(&profile, &codons, 8.0, 1.0, 2.5);
            let nt = align_nt_frame0(&profile, &seq, 8.0, 1.0, 2.5, f64::INFINITY);
            assert!((nt.score - r.score).abs() < 1e-9);
            assert_eq!(nt.sense_span, r.codon_span.map(|(a, b)| (3 * a, 3 * b)));
        }
    }

    #[test]
    fn homopolymer_insertion_gets_the_whole_run_as_its_window() {
        // (MVWFE)x4 + V KKK W + (MVWFE)x4, with one extra A in the AAAAAAAAA run
        let flank = "ATGGTTTGGTTTGAA".repeat(4);
        let mut seq = String::from("CCC");
        seq += &flank;
        seq += "GTT"; // V   sense 63..66
        let run_start = seq.len();
        seq += "AAAAAAAAAA"; // KKK + the extra A
        let run_end = seq.len();
        seq += "TGG"; // W
        seq += &flank;
        let mut aas = b"MVWFE".repeat(4);
        aas.extend_from_slice(b"VKKKW");
        aas.extend(b"MVWFE".repeat(4));
        let profile = single_aa_profile("t", &aas);
        let r = align_nt(&profile, seq.as_bytes(), 8.0, 1.0, 2.5, 10.0);
        assert_eq!(r.frameshifts.len(), 1);
        let f = &r.frameshifts[0];
        assert_eq!(f.n_bases, 1);
        assert!(f.window.0 <= f.sense_pos && f.sense_pos + 1 <= f.window.1);
        // any of the 10 As can be the extra one
        assert_eq!(f.window, (run_start, run_end), "window {:?}", f.window);
    }

    #[test]
    fn reverse_editing_dp_matches_the_reference_implementation() {
        // the score table and the per-codon enumeration must agree with
        // U->C readings in play too, in one pass and in two
        let costs = EditCosts { c_to_u: 2.5, u_to_c: 3.0, u_to_c_stop: 2.0 };
        let mut rng = Lcg(11);
        let mut n_rev = 0;
        for trial in 0..200 {
            let n_codons = 5 + (rng.next() % 40) as usize;
            let (mut seq, profile) = random_gene(&mut rng, n_codons);
            // turn some C into T, as a genome awaiting U->C edits would have
            for i in 0..seq.len() {
                if seq[i] == b'C' && rng.next() % 4 == 0 {
                    seq[i] = b'T';
                }
            }
            let n_junk = (rng.next() % 100) as usize;
            let mut full = random_dna(&mut rng, n_junk);
            full.extend(&seq);
            full.extend(random_dna(&mut rng, n_junk));
            for penalty in [f64::INFINITY, 20.0] {
                for any_frame in [true, false] {
                    let b = align_nt_from_reference(&profile, &full, 8.0, 1.0, costs, penalty, any_frame);
                    for two_pass_min_cells in [usize::MAX, 0] {
                        let a = align_nt_impl(&profile, &full, 8.0, 1.0, costs, penalty, any_frame, two_pass_min_cells);
                        assert_eq!(a.score.to_bits(), b.score.to_bits(), "score, trial {trial}");
                        assert_eq!(format!("{:?}", a.columns), format!("{:?}", b.columns), "columns, trial {trial}");
                    }
                    // and it never scores below C->U editing alone
                    let plain = align_nt_from_reference(&profile, &full, 8.0, 1.0, 2.5, penalty, any_frame);
                    assert!(b.score >= plain.score, "trial {trial}");
                    n_rev += (b.score > plain.score) as usize;
                }
            }
        }
        assert!(n_rev > 50, "reverse editing should matter in this test ({n_rev})");
    }

    #[test]
    fn table_and_packed_pointer_dp_matches_the_reference_implementation() {
        let mut rng = Lcg(7);
        let mut n_shift = 0;
        for trial in 0..400 {
            let n_codons = 5 + (rng.next() % 40) as usize;
            let (mut seq, profile) = random_gene(&mut rng, n_codons);
            // junk either side, the odd N, sometimes an indel
            let n_junk = (rng.next() % 150) as usize;
            let mut full = random_dna(&mut rng, n_junk);
            if trial % 3 == 0 && seq.len() > 12 {
                let at = 3 + (rng.next() as usize) % (seq.len() - 6);
                if trial % 2 == 0 { seq.insert(at, b'A'); } else { seq.remove(at); }
            }
            full.extend(&seq);
            let n_junk = (rng.next() % 150) as usize;
            full.extend(random_dna(&mut rng, n_junk));
            if trial % 5 == 0 {
                let at = (rng.next() as usize) % full.len();
                full[at] = b'N';
            }
            if trial % 7 == 0 {
                let at = (rng.next() as usize) % full.len();
                full[at] = b'c'; // lower case: a non-ACGT byte to the translator
            }
            for penalty in [f64::INFINITY, 20.0, 6.0] {
                for any_frame in [true, false] {
                    let b = align_nt_from_reference(&profile, &full, 8.0, 1.0, 2.5, penalty, any_frame);
                    // one pass (full matrix) and two passes (locate, then a window)
                    for two_pass_min_cells in [usize::MAX, 0] {
                        let a = align_nt_impl(&profile, &full, 8.0, 1.0, 2.5, penalty, any_frame, two_pass_min_cells);
                        assert_eq!(a.score.to_bits(), b.score.to_bits(), "score, trial {trial}");
                        assert_eq!(format!("{:?}", a.columns), format!("{:?}", b.columns), "columns, trial {trial}");
                        assert_eq!(a.sense_span, b.sense_span);
                        assert_eq!(a.frameshifts, b.frameshifts);
                        n_shift += a.frameshifts.len();
                    }
                    if !any_frame {
                        let s = align_nt_frame0_score(&profile, &full, 8.0, 1.0, 2.5, penalty);
                        assert_eq!(s.to_bits(), b.score.to_bits(), "score-only, trial {trial}");
                    }
                }
            }
        }
        assert!(n_shift > 50, "the comparison should exercise frame changes ({n_shift})");
    }
}
