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
    /// Within `ScanParams::end_margin` profile columns of the profile's
    /// start or end: reported, but the gene is not read through it - it
    /// starts after the change, or ends before it in its own frame (the
    /// corrected span, the protein and the edits are those of the gene
    /// without it).
    pub near_end: bool,
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
    /// Literal stop codons inside the gene, each aligned to a profile
    /// column (`*` in `protein`); none unless
    /// `ScanParams::internal_stop_penalty` is finite.
    pub n_internal_stops: usize,
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
    /// `ScanParams::frameshift_penalty` is finite), with those the end
    /// margin kept out of the gene (`FrameshiftCall::near_end`).
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
            internal_stop: self.internal_stop_penalty,
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
    /// and extend the gene to it (0: don't). A profile often lacks a
    /// gene's last residues (boundary-truncated training calls, or a tail
    /// longer in this species than in most), and a gene may end short of
    /// its profile's last columns; either way it ends at its stop.
    pub stop_search: usize,
    /// Where a gene starts: when the first aligned codon is not a start
    /// codon (ATG, or ACG that an edit makes AUG), the nearest ATG within
    /// this many codons is taken instead - upstream, in frame and with no
    /// stop between, or downstream, at a codon the alignment puts in a
    /// column whose references start with Met often enough to score it
    /// positively - and failing an ATG, the nearest ACG the same way
    /// (0: don't). See `move_start`.
    pub start_search: usize,
    /// Cost of a literal stop codon inside a gene, aligned to a profile
    /// column and shown as `*` in the protein
    /// (`codon::EditCosts::internal_stop`); infinite (the default) for none.
    pub internal_stop_penalty: f64,
    /// A frame change or an internal stop with fewer than this many
    /// profile columns before it or after it is not kept (0: kept wherever
    /// it is): the gene starts after it, or ends before it. One near the
    /// end buys the profile's last few columns, and cannot be told from a
    /// gene that simply ends differently there, in its own frame at its
    /// own stop (Arabidopsis atp4, eight columns from the end in a run of
    /// eleven Gs, is annotated to the stop two codons on); one near the
    /// start, from a gene that starts after it. For a short profile, a
    /// quarter of its columns at most.
    pub end_margin: usize,
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
            start_search: 0,
            internal_stop_penalty: f64::INFINITY,
            end_margin: 0,
            reverse_edit_penalty: f64::INFINITY,
            stop_edit_penalty: f64::INFINITY,
            known_edit_penalty: None,
            novel_edit_penalty: None,
        }
    }
}

/// How many times `scan_gene` looks again in a longer window when the
/// gene runs out of the one it has.
const MAX_WIDENINGS: usize = 2;
/// An edit-made stop (CAA, CAG, CGA read as a stop) right after the aligned
/// protein is taken on the strength of the protein ending where the profile
/// does: with more than this many of the profile's last columns unaligned,
/// a literal stop within `stop_search` is looked for first.
const EDITED_STOP_MAX_TRAILING: usize = 2;

/// `profile` aligned in a window, with `params`' frame-change filters; and
/// the frame changes near the profile's ends that the gene is not read
/// through (`drop_breaks_near_ends`).
fn align_window(profile: &Profile, sense_seq: &[u8], params: &ScanParams) -> (align::NtAlignResult, Vec<align::Frameshift>) {
    let mut result =
        align::align_nt(profile, sense_seq, params.gap_open, params.gap_extend, params.edit_costs(), params.frameshift_penalty);
    let filtering = params.frameshift_min_gain > 0.0 || params.frameshift_flank > 0;
    if filtering && !result.frameshifts.is_empty() {
        let single = align::align_nt(profile, sense_seq, params.gap_open, params.gap_extend, params.edit_costs(), f64::INFINITY);
        if result.score - single.score < params.frameshift_min_gain || !frameshifts_supported(profile, &result, params) {
            result = single;
        }
    }
    drop_breaks_near_ends(profile, sense_seq, params, result)
}

/// `result` without the frame changes and internal stops that lie within
/// `params.end_margin` profile columns of either end of the profile, and
/// the frame changes so dropped. If those are all it has, the gene is
/// aligned again in one frame with no stop let through - its own
/// sequence, start to stop. If it has others further in, those are kept,
/// and the columns outside the ones dropped go unaligned.
fn drop_breaks_near_ends(
    profile: &Profile,
    sense_seq: &[u8],
    params: &ScanParams,
    mut result: align::NtAlignResult,
) -> (align::NtAlignResult, Vec<align::Frameshift>) {
    let n = result.columns.len();
    let margin = params.end_margin.min(n / 4);
    let stops: Vec<usize> = result
        .columns
        .iter()
        .enumerate()
        .filter(|(_, c)| matches!(c, NtColumnCall::Aligned { aa: codon::STOP, .. }))
        .map(|(j, _)| j)
        .collect();
    // columns before `head` and from `tail` on are given up
    let (mut head, mut tail, mut kept) = (0, n, 0);
    for c in result.frameshifts.iter().map(|f| f.after_column) {
        if c < margin {
            head = head.max(c);
        } else if c + margin > n {
            tail = tail.min(c);
        } else {
            kept += 1;
        }
    }
    for &j in &stops {
        if j < margin {
            head = head.max(j + 1);
        } else if j + margin >= n {
            tail = tail.min(j);
        } else {
            kept += 1;
        }
    }
    if head == 0 && tail == n {
        return (result, Vec::new());
    }
    let inside = |f: &align::Frameshift| f.after_column >= margin && f.after_column + margin <= n;
    let dropped: Vec<align::Frameshift> = result.frameshifts.iter().filter(|f| !inside(f)).cloned().collect();
    if kept == 0 {
        let costs = codon::EditCosts { internal_stop: f64::INFINITY, ..params.edit_costs() };
        return (align::align_nt(profile, sense_seq, params.gap_open, params.gap_extend, costs, f64::INFINITY), dropped);
    }
    result.frameshifts.retain(inside);
    for (j, column) in result.columns.iter_mut().enumerate() {
        if j < head || j >= tail {
            *column = NtColumnCall::Deleted;
        }
    }
    let ends = |c: &NtColumnCall| match c {
        NtColumnCall::Aligned { sense_start, .. } => Some(*sense_start),
        NtColumnCall::Deleted => None,
    };
    let first = result.columns.iter().find_map(ends);
    let last = result.columns.iter().rev().find_map(ends);
    result.sense_span = first.zip(last).map(|(start, last)| (start, last + 3));
    (result, dropped)
}

/// Bases a window is short by, at its start and at its end (sense), for
/// the profile to fit where the alignment puts it. Two signs of that: the
/// profile columns left unaligned at an end, with the aligned part running
/// to within their length of that end of the window, were never looked
/// for; and counting on from the first aligned column (back from the last)
/// the rest of the profile would run past the window's end (start) - an
/// alignment cut off this way can scatter the remaining columns thinly
/// over what window is left, and leave few of them unaligned at the very
/// end. Zero where the window had room. An anchor is not always the gene:
/// oatkDB's rbcL model is the gene's first 698 bp of 1,430, and 300 bp
/// beyond its hit the window ended at codon 341.
fn window_shortfall(result: &align::NtAlignResult, window_len: usize) -> (usize, usize) {
    let aligned = |(j, c): (usize, &NtColumnCall)| match c {
        NtColumnCall::Aligned { sense_start, .. } => Some((j, *sense_start)),
        NtColumnCall::Deleted => None,
    };
    let n = result.columns.len();
    let first = result.columns.iter().enumerate().find_map(aligned);
    let last = result.columns.iter().enumerate().rev().find_map(aligned);
    let (Some((first_column, start)), Some((last_column, last_start))) = (first, last) else {
        return (0, 0);
    };
    let end = last_start + 3;
    let (leading, trailing) = (first_column, n - 1 - last_column);
    let short = |columns: usize, room: usize| (3 * columns).saturating_sub(room);
    // the profile's last column, counting from the first aligned one; its first, from the last
    let runs_past_end = (start + 3 * (n - first_column)).saturating_sub(window_len);
    let runs_past_start = (3 * (last_column + 1)).saturating_sub(end);
    (short(leading, start).max(runs_past_start), short(trailing, window_len - end).max(runs_past_end))
}

/// Start the gene at a start codon (`ScanParams::start_search`). The
/// alignment starts a gene wherever the profile's first columns match,
/// which is not always a start codon: of 457 single-exon genes in 19
/// RefSeq angiosperm mitogenomes, 42 called genes began at a codon that is
/// none of ATG, ACG, GTG or TTG. Upstream, the gene is lengthened to the
/// start codon (the protein still covers the profile's columns only);
/// downstream, the columns before it are given up.
fn move_start(profile: &Profile, sense_seq: &[u8], params: &ScanParams, result: &mut align::NtAlignResult) {
    let Some((s0, end)) = result.sense_span else {
        return;
    };
    if params.start_search == 0 {
        return;
    }
    let codon_at = |p: usize| (p + 3 <= sense_seq.len()).then(|| [sense_seq[p], sense_seq[p + 1], sense_seq[p + 2]]);
    const STARTS: [[u8; 3]; 2] = [*b"ATG", *b"ACG"];
    if codon_at(s0).is_some_and(|c| STARTS.contains(&c)) {
        return;
    }
    // (codons away, sense position, column to start at; None upstream)
    let mut found: Vec<(usize, usize, Option<usize>, [u8; 3])> = Vec::new();
    for k in 1..=params.start_search {
        let Some(p) = s0.checked_sub(3 * k) else { break };
        let Some(c) = codon_at(p) else { break };
        if codon::translate_codon(c) == codon::STOP || c.iter().any(|b| !b"ACGT".contains(b)) {
            break;
        }
        found.push((k, p, None, c));
    }
    let first_shift = result.frameshifts.iter().map(|f| f.sense_pos).min().unwrap_or(usize::MAX);
    for (j, column) in result.columns.iter().enumerate() {
        let NtColumnCall::Aligned { sense_start: p, .. } = *column else { continue };
        if p <= s0 {
            continue;
        }
        if p >= first_shift || (p - s0) % 3 != 0 || (p - s0) / 3 > params.start_search {
            break;
        }
        if let Some(c) = codon_at(p) {
            if profile.score(j, b'M') > 0.0 {
                found.push(((p - s0) / 3, p, Some(j), c));
            }
        }
    }
    // the nearest ATG (upstream on a tie), else the nearest ACG
    let pick = STARTS.iter().find_map(|start| {
        found.iter().filter(|f| &f.3 == start).min_by_key(|f| (f.0, f.2.is_some())).copied()
    });
    let Some((_, p, column, _)) = pick else {
        return;
    };
    if let Some(j) = column {
        for c in result.columns[..j].iter_mut() {
            *c = NtColumnCall::Deleted;
        }
    }
    result.sense_span = Some((p, end));
}

/// Scan one (species, gene) anchor: extract a flanked genomic window around
/// `anchor_start`/`anchor_end` (0-based half-open, the existing `.ctg.bed`
/// call), try all 3 reading frames against `profile`, and return the
/// best-scoring edit-tolerant alignment. Where the gene runs out of the
/// window (`window_shortfall`), the window is lengthened that way and the
/// gene aligned again.
pub fn scan_gene(
    fasta_idx: &FastaIndex,
    contig: &str,
    anchor_start: i64,
    anchor_end: i64,
    strand: char,
    profile: &Profile,
    params: &ScanParams,
) -> Result<ScanOutcome> {
    // genomic context taken before the anchor's start and after its end
    let (mut before, mut after) = (params.flank, params.flank);
    let mut window = fasta_idx.extract_window(contig, anchor_start - before, anchor_end + after, strand)?;
    let (mut result, mut near_ends) = align_window(profile, &window.2, params);
    for _ in 0..MAX_WIDENINGS {
        let (short_start, short_end) = window_shortfall(&result, window.2.len());
        if short_start == 0 && short_end == 0 {
            break;
        }
        // the gene's start is the window's left end on '+', its right on '-'
        let (more_before, more_after) = if strand == '-' { (short_end, short_start) } else { (short_start, short_end) };
        if more_before > 0 {
            before += more_before as i64 + params.flank;
        }
        if more_after > 0 {
            after += more_after as i64 + params.flank;
        }
        let wider = fasta_idx.extract_window(contig, anchor_start - before, anchor_end + after, strand)?;
        if (wider.0, wider.1) == (window.0, window.1) {
            break; // the contig ends there
        }
        window = wider;
        (result, near_ends) = align_window(profile, &window.2, params);
    }
    let (window_start, window_end, sense_seq) = window;
    move_start(profile, &sense_seq, params, &mut result);
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
    // The gene's stop - see `ScanOutcome::includes_stop`. Right after the
    // aligned protein: a literal one, or, where the protein ends with its
    // profile, one that an edit makes (CAA/CAG/CGA). Failing that, the
    // first literal stop within `stop_search` codons, in frame: the gene
    // runs on past its profile's end, or ends short of the profile's last
    // columns. Failing that too, an edit-made stop right there after all.
    if let Some((_, sense_end)) = result.sense_span {
        let codon_at = |p: usize| (p + 3 <= sense_seq.len()).then(|| [sense_seq[p], sense_seq[p + 1], sense_seq[p + 2]]);
        let literal_stop = |p: usize| codon_at(p).is_some_and(|c| codon::translate_codon(c) == codon::STOP);
        let trailing_deleted = result.columns.iter().rev().take_while(|c| matches!(c, NtColumnCall::Deleted)).count();
        let next = codon_at(sense_end)
            .and_then(|c| codon::codon_hypotheses(c).into_iter().filter(|h| h.aa == codon::STOP).min_by_key(|h| h.n_edits));
        let further_on = (1..=params.stop_search).map(|i| sense_end + 3 * i).find(|&p| literal_stop(p));
        // (sense position the gene ends at, the edit that makes its stop)
        let stop = match next {
            Some(hyp) if hyp.n_edits == 0 => Some((sense_end + 3, None)),
            Some(hyp) if trailing_deleted <= EDITED_STOP_MAX_TRAILING || further_on.is_none() => {
                // A stop that needs an edit, with a literal one straight
                // after it: the codon is the protein's last residue
                // (Gln/Arg, missing from the profile) and the literal stop
                // the gene's.
                let p = sense_end + 3;
                if params.stop_search > 0 && literal_stop(p) {
                    Some((p + 3, None))
                } else {
                    Some((p, Some(hyp)))
                }
            }
            _ => further_on.map(|p| (p + 3, None)),
        };
        if let Some((stop_end, edited)) = stop {
            let (ext_start, ext_end) = sense_span_to_genomic(window_start, window_end, strand, sense_end, stop_end);
            if strand == '+' {
                corrected_end = ext_end;
            } else {
                corrected_start = ext_start;
            }
            includes_stop = true;
            if let Some(hyp) = edited {
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
    let mut frameshifts: Vec<FrameshiftCall> = result
        .frameshifts
        .iter()
        .map(|f| (f, false))
        .chain(near_ends.iter().map(|f| (f, true)))
        .map(|(f, near_end)| {
            let (start, end) =
                sense_span_to_genomic(window_start, window_end, strand, f.sense_pos, f.sense_pos + f.n_bases);
            let (ws, we) = sense_span_to_genomic(window_start, window_end, strand, f.window.0, f.window.1);
            FrameshiftCall { start, end, n_bases: f.n_bases, profile_column: f.after_column, window_start: ws, window_end: we, near_end }
        })
        .collect();
    frameshifts.sort_by_key(|f| f.start);
    edits.sort_by_key(|e| e.genomic_pos);

    let n_internal_stops = protein.bytes().filter(|&b| b == codon::STOP).count();
    Ok(ScanOutcome {
        frame,
        score: result.score,
        corrected_start,
        corrected_end,
        protein,
        edits,
        n_columns_aligned: n_aligned,
        n_columns_total: profile.len(),
        n_internal_stops,
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

    fn profile_of(aas: &[u8]) -> Profile {
        let columns = aas
            .iter()
            .map(|target| {
                let mut col = [-5.0f64; 20];
                col[AA_ORDER.iter().position(|x| x == target).unwrap()] = 5.0;
                col
            })
            .collect();
        Profile { gene: "test".to_string(), columns, edited: None }
    }

    fn revcomp(seq: &str) -> String {
        seq.bytes().rev().map(|b| match b { b'A' => 'T', b'C' => 'G', b'G' => 'C', _ => 'A' }).collect()
    }

    #[test]
    fn a_gene_that_runs_out_of_its_window_is_looked_for_in_a_longer_one() {
        // the anchor is the gene's first three codons of nine (a model of
        // part of the gene), and with a flank of 3 the window ends a codon
        // beyond it
        let gene = "ATGAAAGTTTGGTTTGAATATCATCCG"; // M K V W F E Y H P
        let profile = profile_of(b"MKVWFEYHP");
        let contig = format!("CCC{gene}TAACCCCCCCCCCCC");
        let params = ScanParams { flank: 3, ..ScanParams::default() };
        let out = scan_gene(&fasta_with("ctg", &contig), "ctg", 3, 12, '+', &profile, &params).unwrap();
        assert_eq!((out.n_columns_aligned, out.corrected_start, out.corrected_end, out.includes_stop), (9, 3, 33, true));
        assert_eq!(out.protein, "MKVWFEYHP");
        // the same gene on the other strand: its start is the anchor's right end
        let n = contig.len() as i64;
        let out = scan_gene(&fasta_with("ctg", &revcomp(&contig)), "ctg", n - 12, n - 3, '-', &profile, &params).unwrap();
        assert_eq!((out.n_columns_aligned, out.corrected_start, out.corrected_end, out.includes_stop), (9, n as usize - 33, n as usize - 3, true));
        // an anchor on the gene's last three codons: the window grows the other way
        let out = scan_gene(&fasta_with("ctg", &contig), "ctg", 21, 30, '+', &profile, &params).unwrap();
        assert_eq!((out.n_columns_aligned, out.corrected_start, out.corrected_end), (9, 3, 33));
        // with two codons more in it than the profile has, the gene takes a second look
        let long = format!("CCC{}GGTGCA{}TAACCCCCCCCCCCC", &gene[..12], &gene[12..]);
        let out = scan_gene(&fasta_with("ctg", &long), "ctg", 3, 12, '+', &profile, &params).unwrap();
        assert_eq!((out.n_columns_aligned, out.corrected_start, out.corrected_end, out.includes_stop), (9, 3, 39, true));
        // a gene that really is its first three codons, with room to see so, is left alone
        let short = format!("CCC{}TAACCCCCCCCCCCCCCCCCCCCCCCCCCC", &gene[..9]);
        let out = scan_gene(&fasta_with("ctg", &short), "ctg", 3, 12, '+', &profile, &ScanParams { flank: 30, ..params }).unwrap();
        assert_eq!((out.n_columns_aligned, out.corrected_start, out.corrected_end, out.includes_stop), (3, 3, 15, true));
        // and so is one the contig ends in
        let cut = format!("CCC{}", &gene[..12]);
        let out = scan_gene(&fasta_with("ctg", &cut), "ctg", 3, 12, '+', &profile, &params).unwrap();
        assert_eq!((out.n_columns_aligned, out.corrected_end), (4, 15));
    }

    #[test]
    fn a_gene_ending_short_of_its_profile_ends_at_its_stop() {
        // M K V W F E, a seventh residue the profile does not have, and
        // the stop: the profile's last three columns stay unaligned
        let profile = profile_of(b"MKVWFEYYY");
        let idx = fasta_with("ctg", "ATGAAAGTTTGGTTTGAACCCTAAGGGGGGGGGGGG");
        let params = ScanParams { flank: 30, stop_search: 2, ..ScanParams::default() };
        let out = scan_gene(&idx, "ctg", 0, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.n_columns_aligned, out.includes_stop, out.corrected_end), (6, true, 24));
        assert_eq!(out.protein, "MKVWFE---");
        // CAA there could be a stop made by an edit, but with the profile's
        // end unaligned the literal stop two codons on is the gene's
        let idx = fasta_with("ctg", "ATGAAAGTTTGGTTTGAACAAGGGTAAGGGGGGGGG");
        let out = scan_gene(&idx, "ctg", 0, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end, out.edits.len()), (true, 27, 0));
        // with none in reach, the edit-made one it is
        let out = scan_gene(&idx, "ctg", 0, 18, '+', &profile, &ScanParams { stop_search: 1, ..params }).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end, out.edits.len()), (true, 21, 1));
        // where the protein ends with its profile, an edit-made stop stands
        let profile = profile_of(b"MKVWFE");
        let out = scan_gene(&idx, "ctg", 0, 18, '+', &profile, &params).unwrap();
        assert_eq!((out.includes_stop, out.corrected_end, out.edits.len()), (true, 21, 1));
    }

    #[test]
    fn frameshift_end_margin_keeps_a_gene_in_its_own_frame_at_its_end() {
        // (MVWFE)x4 W, an extra base, H D Y R N and a stop in that other
        // frame; in the gene's own frame the stop is one codon after the W
        let seq = format!("{}TGGCCATGATTATCGTAATTAGCCCCCC", "ATGGTTTGGTTTGAA".repeat(4));
        let mut aas = b"MVWFE".repeat(4);
        aas.extend_from_slice(b"WHDYRN");
        let profile = profile_of(&aas);
        let idx = fasta_with("ctg", &seq);
        let params = ScanParams { flank: 0, frameshift_penalty: 10.0, stop_search: 10, ..ScanParams::default() };
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &params).unwrap();
        assert_eq!((out.frameshifts.len(), out.corrected_end), (1, 82));
        // five columns follow the frame change: kept with a margin of 5, not of 6
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { end_margin: 5, ..params }).unwrap();
        assert_eq!((out.frameshifts.len(), out.corrected_end), (1, 82));
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { end_margin: 6, ..params }).unwrap();
        // it is still reported, at the same place, as one the gene is not read through
        assert_eq!(out.frameshifts.iter().map(|f| (f.start, f.near_end)).collect::<Vec<_>>(), [(63, true)]);
        assert_eq!((out.includes_stop, out.corrected_end), (true, 69), "its own frame, to the stop there");
        assert_eq!(out.protein, format!("{}W-----", "MVWFE".repeat(4)));
        // with a frame change further in as well, that one is kept and the gene's end left unaligned
        let seq = format!("{}C{}", "ATGGTTTGGTTTGAA".repeat(4), seq);
        let mut aas = b"MVWFE".repeat(8);
        aas.extend_from_slice(b"WHDYRN");
        let profile = profile_of(&aas);
        let idx = fasta_with("ctg", &seq);
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &params).unwrap();
        assert_eq!(out.frameshifts.len(), 2);
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { end_margin: 6, ..params }).unwrap();
        assert_eq!(out.frameshifts.iter().map(|f| (f.profile_column, f.near_end)).collect::<Vec<_>>(), [(20, false), (41, true)]);
        assert_eq!((out.includes_stop, out.corrected_end), (true, 130));
        assert_eq!(out.protein, format!("{}W-----", "MVWFE".repeat(8)));
    }

    #[test]
    fn a_stop_with_the_gene_going_on_beyond_it_is_an_internal_stop() {
        // (MVWFE)x3, a stop where the profile has K, (MVWFE)x3 and the gene's stop
        let unit = "ATGGTTTGGTTTGAA";
        let mut aas = b"MVWFE".repeat(3);
        aas.push(b'K');
        aas.extend(b"MVWFE".repeat(3));
        let profile = profile_of(&aas);
        let seq = format!("CCC{}TGA{}TAGCCC", unit.repeat(3), unit.repeat(3));
        let idx = fasta_with("ctg", &seq);
        let params = ScanParams { flank: 0, stop_search: 10, ..ScanParams::default() };
        // no stop let through: one side of it is the gene
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &params).unwrap();
        assert_eq!((out.n_internal_stops, out.n_columns_aligned), (0, 15));
        // let through at a price, it shows, and the gene is whole
        let params = ScanParams { internal_stop_penalty: 20.0, ..params };
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &params).unwrap();
        assert_eq!((out.n_internal_stops, out.corrected_start, out.corrected_end, out.includes_stop), (1, 3, 99, true));
        assert_eq!(out.protein, format!("{0}*{0}", "MVWFE".repeat(3)));
        // fifteen columns either side of it: kept with that margin, and with a larger one
        // capped at a quarter of this short profile; a frame change would be too
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { end_margin: 30, ..params }).unwrap();
        assert_eq!(out.n_internal_stops, 1);
        // a stop five columns from the profile's end is the gene's own stop
        let mut aas = b"MVWFE".repeat(6);
        aas.push(b'K');
        aas.extend(b"MVWFE");
        let profile = profile_of(&aas);
        let seq = format!("CCC{}TGA{}TAGCCC", unit.repeat(6), unit);
        let idx = fasta_with("ctg", &seq);
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &params).unwrap();
        assert_eq!((out.n_internal_stops, out.corrected_end), (1, 114));
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { end_margin: 6, ..params }).unwrap();
        assert_eq!((out.n_internal_stops, out.corrected_start, out.corrected_end, out.includes_stop), (0, 3, 96, true));
        // and one five columns from its start is before the gene
        let mut aas = b"MVWFE".to_vec();
        aas.push(b'K');
        aas.extend(b"MVWFE".repeat(6));
        let profile = profile_of(&aas);
        let seq = format!("CCC{}TGA{}TAGCCC", unit, unit.repeat(6));
        let idx = fasta_with("ctg", &seq);
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { end_margin: 6, ..params }).unwrap();
        assert_eq!((out.n_internal_stops, out.corrected_start, out.corrected_end), (0, 21, 114));
    }

    #[test]
    fn a_gene_starts_at_a_start_codon() {
        // the profile is K V W F E Y; the gene starts at the ATG two codons
        // before the K (GGC between), and the alignment would start at the K
        let profile = profile_of(b"KVWFEY");
        let seq = "CCCTAACCCATGGGCAAAGTTTGGTTTGAATATTAG";
        let params = ScanParams { flank: 0, stop_search: 10, ..ScanParams::default() };
        let idx = fasta_with("ctg", seq);
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &params).unwrap();
        assert_eq!((out.corrected_start, out.corrected_end), (15, 36));
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { start_search: 2, ..params }).unwrap();
        assert_eq!((out.corrected_start, out.corrected_end, out.n_columns_aligned), (9, 36, 6));
        assert_eq!(out.protein, "KVWFEY");
        // one codon is not far enough
        let out = scan_gene(&idx, "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { start_search: 1, ..params }).unwrap();
        assert_eq!(out.corrected_start, 15);
        // nor past a stop
        let seq = "CCCCCCCCCATGTAAAAAGTTTGGTTTGAATATTAG";
        let out = scan_gene(&fasta_with("ctg", seq), "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { start_search: 5, ..params }).unwrap();
        assert_eq!(out.corrected_start, 15);
        // downstream: an ATG at a column whose references have Met; the columns before it are given up
        let profile = profile_of(b"KMVWFEY");
        let seq = "CCCCCCCCCAAAATGGTTTGGTTTGAATATTAG";
        let out = scan_gene(&fasta_with("ctg", seq), "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { start_search: 5, ..params }).unwrap();
        assert_eq!((out.corrected_start, out.corrected_end, out.n_columns_aligned), (12, 33, 6));
        assert_eq!(out.protein, "-MVWFEY");
        // not at a column the references never start at
        let profile = profile_of(b"KKVWFEY");
        let seq = "CCCCCCCCCAAAAAAGTTTGGTTTGAATATTAG";
        let out = scan_gene(&fasta_with("ctg", seq), "ctg", 0, seq.len() as i64, '+', &profile, &ScanParams { start_search: 5, ..params }).unwrap();
        assert_eq!(out.corrected_start, 9);
        // an ACG (an edit away from AUG) counts as a start already
        let profile = profile_of(b"KVWFEY");
        let seq = "CCCCCCATGACGAAAGTTTGGTTTGAATATTAG";
        let profile2 = profile_of(b"TKVWFEY");
        let out = scan_gene(&fasta_with("ctg", seq), "ctg", 0, seq.len() as i64, '+', &profile2, &ScanParams { start_search: 5, ..params }).unwrap();
        assert_eq!(out.corrected_start, 9);
        let _ = profile;
        // and on the minus strand, upstream is to the right
        let seq = "CCCTAACCCATGGGCAAAGTTTGGTTTGAATATTAG";
        let rc = revcomp(seq);
        let n = rc.len();
        let profile = profile_of(b"KVWFEY");
        let out = scan_gene(&fasta_with("ctg", &rc), "ctg", 0, n as i64, '-', &profile, &ScanParams { start_search: 2, ..params }).unwrap();
        assert_eq!((out.corrected_start, out.corrected_end), (0, n - 9));
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
