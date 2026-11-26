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
   A profile often lacks a gene's last residues (it is built from gene
   calls, some cut short; or the gene's tail is longer in this species
   than in most), and a gene may end short of its profile's last columns.
   So when no stop follows directly, `scan-batch` looks up to
   `--stop-search` codons (default 50; 0 for the adjacent codon only)
   further on, in frame, for a literal stop and ends the gene there; and
   a codon that would be a stop only after an edit (`CAA`/`CAG`/`CGA`)
   with a literal stop straight after it is taken as the last residue.
   The protein column still covers the profile's columns only. See "Where
   a gene ends" below.
5. The anchor is not always the whole gene (oatkDB's rbcL model is the
   gene's first 698 bp of 1,430, and 300 bp beyond its hit the window
   ended at codon 341). Where the alignment runs into an end of its
   window with profile columns left over, or counting on from its first
   aligned column the profile would run past the window, the window is
   lengthened that way and the gene aligned again (twice at most).

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

## Where a gene starts

The alignment starts a gene wherever the profile's first columns match,
which is not always at a start codon: of 457 single-exon genes in 19
RefSeq angiosperm mitogenomes, 42 called genes began at a codon that is
none of ATG, ACG, GTG or TTG. With `--start-search N` (default 10; 0 to
start where the alignment does), a gene whose first aligned codon is
neither ATG nor ACG (which an edit makes AUG) starts at the nearest ATG
within N codons instead: upstream, in frame and with no stop between
(the gene is lengthened; the protein still covers the profile's columns
only), or downstream, at a codon the alignment puts in a column where
the references' Met scores positively (the columns before it are given
up). Failing an ATG, the nearest ACG is taken the same way; failing
both, the start is left alone.

On the single-exon genes of 19 RefSeq angiosperm mitogenomes and 13
RefSeq plastomes, starts within 3 bp of the annotated ones went from 384
to 401 of 457 (mito) and from 528 to 574 of 626 (plastid; 564 to 613 of
671 counting the Pinus plastome, which the harness had left out), with one
start that was right made wrong (Arabidopsis ccmFn2, annotated at a GTG
six codons from an ATG). At 5, 15, 20 and 30 codons: 393/559, 407/576,
409/577 and 409/576, with 1, 6, 8 and 9 made wrong - genes annotated
with another start codon (TTG, GTG, ATT, ATA) and an ATG near them.
What is left is mostly genes whose annotated start lies far from the
conserved part (atp6 has N-terminal extensions of 18 to 146 codons in
these genomes) or differs between records (rps4, mttB, rpl16).

`--upstream-start GENES` (0.2.1; comma-separated, matched to the
manifest's gene name before any `|` or to the profile's file name) is
for such a gene: one whose N-terminal extension is too variable to align
and is annotated from the first ATG of its open reading frame. It starts
at the most upstream ATG of the frame, up to 300 codons before the
aligned start with no stop between, unless the gene starts at an ACG (a
start an edit makes); with no ATG upstream, `--start-search` applies as
usual. The frame is followed in the contig after the gene is aligned, so
the alignment, its score, edits and frame changes, and the gene's end
are as without the option; only its start moves. With `--upstream-start
atp6` on the same 19 mitogenomes, atp6 starts within 3 bp went from 8 to
14 of 19 and mito starts from 401 to 407 of 457, none made wrong; of the
five atp6 still wrong, four now start further upstream than RefSeq does,
and one at an ACG that RefSeq starts after. Over 1,785 atp6 copies in 836
DToL angiosperm mitogenomes, 861 start further upstream (median 87
codons, 90% within 179), and nothing else about them changes.

## Where a gene ends

A reading frame ends at its first stop codon. The aligner used to treat
a stop as one more codon it could leave out, so wherever a gene ended a
few residues before its profile did, the alignment stepped over the stop
to match the profile's last columns (maize cob ends `...KYYTE*`, and was
read on through two stops), and likewise before a gene's start. A stop
codon now cannot lie in an insertion (unless stops are read through,
`--stop-edit-penalty`). Three rules follow from that:

- **The stop.** A literal stop right after the aligned protein ends the
  gene; so does one an edit makes (`CAA`/`CAG`/`CGA`) there, provided
  the protein ends with its profile (at most two of the profile's last
  columns unaligned). Otherwise the first literal stop within
  `--stop-search` codons does, and only failing that an edit-made one.
- **Internal stops** (`--internal-stop-penalty`, default 20; `inf` for
  none). A gene whose protein goes on matching beyond a stop - a
  pseudogene, a stop a lineage edits away, an error - has the stop
  aligned to a profile column at this cost: `*` in the protein, counted
  in `n_internal_stops`. Without it such a gene is cut at the stop, or
  the alignment gets round the stop with two frame changes; the cost
  should not exceed that of one frame change.
- **The margin** (`--end-margin`, default 25 profile columns; a quarter
  of a short profile at most). The gene is not read through an internal
  stop or a frame change with fewer columns than this before it or
  after it: it starts after it, or ends before it in its own frame at
  its own stop. One that buys only the profile's last few columns
  cannot be told from a gene that ends differently there (Arabidopsis
  atp4 changes frame eight columns from the end, in a run of eleven Gs,
  to match the other angiosperms' tail; it is annotated to the stop two
  codons on). If the gene has no other breaks it is aligned again in one
  frame; if it has, those are kept and the columns outside go unaligned.
  The frame change itself is still reported, in `<frameshifts-out>` with
  `near_end` true (and not counted in `n_frameshifts`): an indel there is
  as worth checking against the reads as any other. Dropping them
  outright cost 12 points of the 1-bp indels planted in real genes that
  are found (95% to 83% of 1,500 at a margin of 25; 94%, 90% and 86% at
  10, 15 and 20).

On the single-exon genes of 19 RefSeq angiosperm mitogenomes and 13
RefSeq plastomes (`poa`'s `build/rnaedit/ends`), with the same profiles
before and after: genes ending within 3 bp of the annotated stop went
from 397 to 443 of 457 (mito) and from 558 to 621 of 626 (plastid), and
rbcL from 341 of 474 codons to all of them; one end that was right was
made wrong (an edit-made stop a codon beyond the profile's end, where a
literal one 13 codons on is now taken). Frame changes called in those
mito genes went from 23 to 3. Margins of 15, 20, 25 and 30 gave 436,
440, 443 and 443 mito stops; internal-stop costs of 15 to 30 the same
ends; stop searches of 10, 20, 30, 50 and 100 codons 583, 596, 611, 617
and 618 plastid stops (before the longer windows). A gene that ends
before its profile's tail no longer counts as aligned over the whole
profile, which it used to by the step over the stop.

## Reverse (U-to-C) editing

Seed plants, mosses and liverworts edit C to U only. Ferns, quillworts
and hornworts also edit U to C, and heavily: a third of the edits in a
fern mitochondrion, and nearly half of those turn a genomic stop codon
into Gln or Arg, about nine stops per gene in Polypodiales. With C-to-U
alone each such stop is a column the alignment has to skip, or looks like
the end of the gene.

Two options, both off (`inf`) by default:

- `--stop-edit-penalty`: reading an in-frame TAA or TAG as Gln, or TGA as
  Arg.
- `--reverse-edit-penalty`: any other genomic T read as C.

They are separate because a stop inside a gene is strong evidence of an
edit and an ordinary T is not. Measured on 702 single-exon mito gene
copies from 30 DToL fern assemblies, against 4,380 copies from 60
angiosperm assemblies in which any U-to-C call is false:

- `--stop-edit-penalty 2.5` (the cost of a C-to-U edit): fern copies
  aligned over 99% of their profile go from 152 to 292, frameshift calls
  from 277 to 113, and U-to-C comes to 32% of mito and 6% of plastid
  edits - 33-35% and 1-10% are measured in fern transcriptomes. At 4 and
  6 fewer copies align (289, 285).
- Adding `--reverse-edit-penalty` 8 or 12 changes nothing for the ferns.
  One flat cost of 2.5 for every reverse edit gives the same 292 but
  calls U-to-C for 55% of fern mito edits and puts two false ones in
  every angiosperm gene copy.
- In the angiosperm copies, stop read-through at 2.5 changes 461
  proteins: stops in truncated or pseudogene copies read through. It is
  for lineages that edit this way, not a safe default.

With either option off the output is what it was (every table for 90,102
DToL gene copies unchanged). The edits table has an `edit` column, `C>U`
or `U>C`.

## Editing sites

A profile built from RNA-edited proteins can say where they are edited: a
21st column, `edit`, holding for each profile column the share of its
reference proteins whose residue there an edit made (0: none). Without
it, a C-to-U edit costs `--edit-penalty` at every column, and wherever a
species has a residue of its own that a C-to-U would turn into the
profile's, an edit is called: in fast-evolving genes (plastid rpoC2,
matK) most calls were of that kind. With it, two more costs apply to
that profile, each `--edit-penalty` unless given:

- `--known-edit-penalty`: a C-to-U edit at a column edited in the
  references.
- `--novel-edit-penalty`: one at any other column.

Measured on 107 DToL angiosperms with RNA-seq, each annotated with
profiles built without its own genus: at 2.5 everywhere, 40% of the
edits called in single-exon plastid genes and 93% of those in
mitochondrial ones were in the RNA. With 10 at the other columns and 1.5
at the known ones, 94% and 97% are, with more correct plastid calls
than before (1,494 against 1,444) and 9% more mitochondrial ones, which
find 77% of the C-to-U sites the reads show, from 71%.

The profile's builder decides what counts as a site; poa asks for the
edit in two genera, since RNA-seq from another plant than the assembly
shows the differences between the two as edits. Stop codons made by an
edit are not columns and keep `--edit-penalty`. A profile without the
column is read as before, and the output is unchanged for it.

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
frame changes in the gene, `n_internal_stops` (the last column) the stop
codons aligned inside it (`*` in the protein), and `frame` is the frame
(0-2, within the window) of its first aligned codon. `<edits-out>`: one row per
predicted edit site (genomic position, which profile column - or the
profile length itself, for a predicted edit-created stop - and the
resulting amino acid). `<frameshifts-out>` (optional): one row per
frame change - `start`/`end` (0-based half-open genomic span of the
skipped bases), `n_bases` (1 or 2), `profile_column` (the column after
it), `window_start`/`window_end` (the window, containing start/end) and
`near_end` (true for one within `--end-margin` of the profile's start
or end, which the gene's span and protein do not follow).

## Profile format

Plain TSV: a header row of the 20 amino acid one-letter codes (any column
order, one column per code, exactly once each), then one row per profile
(alignment) column of per-amino-acid scores (log-odds or any consistent
scoring scheme - `orfedit` just sums them, it doesn't care how they were
derived). See `plant_organellar_database`'s
`analysis/editing/src/00_build_reference_profiles.py` for how these are
built from that dataset's own annotated orthologs. A 21st column headed
`edit` is optional (see Editing sites): a number per profile column, above
0 where the reference proteins are edited.

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
`--frameshift-penalty`, `--stop-search`, `--internal-stop-penalty`, `--end-margin`, `--stop-edit-penalty`, `--reverse-edit-penalty`, and `--flank` (bp of genomic context to search around each anchor, default 300; lengthened for a gene that runs out of it).

## Scope

Single-exon genes only. Genes that are trans-spliced across separate
transcript fragments (common for plant mitochondrial `nad1/2/5`, `rps3`)
need exon reconstruction first, which is a separate, not-yet-built tool -
`orfedit`'s alignment core is designed to be reusable there too (validating
a candidate joined-exon ORF is the same edit-tolerant profile-alignment
problem, just with the genome-side input being a concatenation of exon
fragments instead of one contiguous window).
