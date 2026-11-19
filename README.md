# orfedit

A homology-guided, RNA-editing-tolerant gene boundary corrector for plant
organellar (mitochondrial/plastid) genomes.

## The problem

Plant mitochondrial (and, more sparingly, plastid) mRNAs undergo extensive
C-to-U RNA editing after transcription - including, routinely, at start and
stop codons (`ACG` -> `AUG` start creation; a genomic sense codon like
`CGA`/`CAA`/`CAG` that only becomes a stop, `UGA`/`UAA`/`UAG`, after
editing). A DNA-level ORF finder or HMM gene caller has no way to see this:
it either reads straight through an edited stop (over-extending into
downstream sequence) or misses an edited start (truncating the N-terminus,
or missing the gene call outright).

`orfedit` doesn't predict editing sites from first principles (that's what
tools like PREP-Mt do, from a curated, trained multi-species model). Instead
it takes a homology profile built from a large set of already-annotated
orthologs of the same gene, and asks: for this one genomic locus, which
combination of "read this C literally" vs. "read this C as if edited to U"
best reconstructs a full-length protein matching that profile? It's an
edit-tolerant profile-to-DNA aligner, not an editing-site predictor per se -
the predicted edits are a side effect of finding the best-scoring alignment,
not the primary target.

## How it works

For one (genomic locus, reference profile) pair:

1. Extract a flanked genomic window around the existing (approximate,
   DNA-level) gene call.
2. Run a semi-global affine-gap DP alignment of the amino acid profile
   against the window's codons in all 3 reading frames at once (the DP runs
   over sense bases; a codon is the 3 bases ending at a position) - free
   end-gaps on the genome side (so the true start/stop can be anywhere in
   the window), full-length required on the profile side (so a real,
   complete protein is what's being searched for). Between two codons the
   alignment may also skip 1 or 2 bases, at a cost
   (`--frameshift-penalty`), moving the rest of the gene into another
   frame - see Frameshifts below.
3. At every codon, the substitution score used is the *best* of: translate
   it literally, or translate it with any subset of its C positions
   hypothesized as edited to T, each hypothesized edit costing a small fixed
   penalty (`--edit-penalty`) so the aligner doesn't invent edits gratuitously.
4. Traceback gives the corrected start/stop coordinates, the edit-corrected
   protein, and the list of individual predicted edit positions. One extra
   step after traceback: profiles are built from stop-stripped protein
   translations, so the aligned span alone is always one codon short of a
   stop-inclusive gene-boundary convention (e.g. `.ctg.bed`-style calls);
   `orfedit` checks the codon immediately after the aligned protein and, if
   it's a genuine (possibly edit-created) stop, folds it into the corrected
   span (`ScanOutcome::includes_stop`) rather than reporting a spurious
   "correction" on every single call regardless of real editing.

Interior profile columns can only ever resolve to a real amino acid or a
deletion, never to a stop/unrecognised translation: `Profile::score` scores
those as `NEG_INFINITY`, not a finite penalty. This was a real bug once
(a finite -10.0 penalty let a stop-codon reading win an interior column
when every real-residue hypothesis scored even worse, leaking literal `*`
characters into output protein sequences), caught downstream by
`transsplice`'s validation, fixed here at the root rather than worked
around in a caller - see `src/profile.rs` and the regression test
`align::tests::unrescuable_interior_stop_is_deleted_not_aligned`.

See `src/align.rs` for the DP itself and `src/codon.rs` for the edit
hypothesis enumeration - both have unit tests covering the two canonical
cases (`ACG`->`ATG` start rescue, an edited internal stop) plus a
stability/specificity check (no editing needed -> no change reported).
`src/lib.rs`'s tests cover the trailing-stop folding above.

## Frameshifts

A frameshifting indel - an assembly error such as a HiFi homopolymer
miscount, or a real pseudogene - leaves a single-frame alignment able to
match only one side of it, the rest of the profile deleted. With frame
changes allowed, the alignment carries on in the new frame and reports
where it changed: the 1 or 2 skipped bases. Skipping 1 base means the
gene has one base more than its reading frame wants (a 1-bp insertion, or
a 2-bp deletion); skipping 2 means one base fewer (a 1-bp deletion, or a
2-bp insertion).

Where the change could equally well sit is often wider than the skipped
bases - anywhere in a homopolymer run, or along a stretch both frames
read about as well - so each frame change also gets a **window**: every
breakpoint within the ungapped 30 columns either side whose score comes
within 2 of the best one. The extra A in `GTT AAAAAAAAAA TGG` gets the
whole run.

`--frameshift-penalty` (default 20) is the cost of one frame change;
`inf` turns them off and gives exactly the 0.1 result. Two optional
filters revert a gene to its single-frame alignment: `--frameshift-flank`
(> 0: that many soundly aligned profile columns needed each side, with
`--frameshift-max-deleted` / `--frameshift-min-flank-score`) and
`--frameshift-min-gain` (the frame changes must beat the best single
frame by that much). Both are off: they mostly remove calls in partial
gene copies and ones that recur across assemblies, and cost recall.

Validation on the DToL plant organelle dataset (90,375 gene copies from
2,122 assemblies; `poa`'s `build/fs_val` scripts):

- Frame changes off: output identical to 0.1.0. On: every gene without a
  frame change is aligned exactly as with them off.
- 1,500 1-bp indels planted at random in the coding sequence of real gene
  copies are found (at the indel, or with it inside the window) in 92%:
  96% in plastid genes, 88% in mito ones, 94% in genes' best copies
  aligning 90% or more of the profile; 99.4% with the right type. The
  misses are mostly short genes (rps19, atp9, rps14) and indels near a
  gene end, where too little profile is left on one side to pay for a
  frame change.
- Across the dataset, a mean of 2 calls per assembly (median 1); 79 calls
  in a gene's best copy at a place no other assembly has one - the
  candidates for assembly errors.
- Penalty 30 instead: 90% found, a mean of 1.1 calls per assembly.
  The 10/2/1.5 flank filter at penalty 30: 81% found.

Speed: the alignment looks codon scores up in a table built once per
profile (125 codon classes per column), keeps one byte of traceback per
cell, and for a large search (8 million cells or more - a whole contig,
as `transsplice` asks for) first finds where the best alignment lies
without any traceback, then aligns only that stretch. Results are the
same, bit for bit (a unit test compares against the original
implementation; all three output tables for the 90,375 DToL copies are
unchanged), in 55 s instead of 6 min 56 s at 16 threads.
`align::align_nt_frame0_score` gives an alignment's score alone.

## Output (`scan-batch`)

`<gene-calls-out>`: one row per manifest row - `raw_start`/`raw_end` (echoed
from the manifest) alongside `corrected_start`/`corrected_end`,
`includes_stop`, `n_edits`, alignment `score`, and the edit-corrected
protein (`-` = a profile column with no aligned genomic codon). `changed`
means `n_edits > 0` ("editing was invoked"), deliberately not "coordinates
differ at all" - a profile-alignment boundary essentially never lands on
the exact same bp as an independently-called boundary even with zero real
editing, so raw coordinate (in)equality is not a usable "did anything
biological happen" signal; compare `raw_*`/`corrected_*` yourself if you
want boundary drift regardless of editing. `n_frameshifts` counts the
frame changes in the gene, and `frame` is the frame (0-2, within the
window) of its first aligned codon. `<edits-out>`: one row per
predicted edit site (genomic position, which profile column - or the
profile length itself, for a predicted edit-created stop - and the
resulting amino acid). `<frameshifts-out>` (optional): one row per
frame change - `start`/`end` (0-based half-open genomic span of the
skipped bases), `n_bases` (1 or 2), `profile_column` (the column after
it) and `window_start`/`window_end` (the window, containing start/end).

## Profile format

Plain TSV: a header row of the 20 amino acid one-letter codes (any column
order, one column per code, exactly once each), then one row per profile
(alignment) column of per-amino-acid scores (log-odds or any consistent
scoring scheme - `orfedit` just sums them, it doesn't care how they were
derived). See `plant_organellar_database`'s
`analysis/editing/src/00_build_reference_profiles.py` for how these are
built from that dataset's own annotated orthologs.

## Usage

```
orfedit scan-batch \
  --manifest manifest.tsv \
  --gene-calls-out edited_gene_calls.tsv \
  --edits-out predicted_edits.tsv \
  --frameshifts-out frameshifts.tsv \
  --threads 16
```

`manifest.tsv` header (tab-separated):
`species organelle gene fasta profile contig start end strand`
(`start`/`end` are 0-based half-open genomic coordinates of the existing
DNA-level gene call - `orfedit` extends/corrects outward from there, it
doesn't discover genes de novo.)

Tunable scoring: `--gap-open`, `--gap-extend`, `--edit-penalty`,
`--frameshift-penalty`, and `--flank` (bp of genomic context to search around each anchor, default 300).

## Scope

Single-exon genes only. Genes that are trans-spliced across separate
transcript fragments (common for plant mitochondrial `nad1/2/5`, `rps3`)
need exon reconstruction first, which is a separate, not-yet-built tool -
`orfedit`'s alignment core is designed to be reusable there too (validating
a candidate joined-exon ORF is the same edit-tolerant profile-alignment
problem, just with the genome-side input being a concatenation of exon
fragments instead of one contiguous window).
