//! `orfedit` - homology-guided, C->U-edit-tolerant gene boundary correction
//! for plant organellar genomes. See README.md for the algorithm and the
//! `plant_organellar_database` repo's `analysis/editing/` module for how
//! this is invoked in practice.

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use orfedit::fasta::FastaIndex;
use orfedit::profile::Profile;
use orfedit::{scan_gene, ScanParams};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "orfedit", version, about)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Scan every (species, gene) row in a manifest TSV in parallel.
    ScanBatch {
        /// TSV with header: species organelle gene fasta profile contig start end strand
        #[arg(long)]
        manifest: PathBuf,
        #[arg(long)]
        gene_calls_out: PathBuf,
        #[arg(long)]
        edits_out: PathBuf,
        /// Also write one row per frameshift found in a corrected gene
        #[arg(long)]
        frameshifts_out: Option<PathBuf>,
        #[arg(long, default_value_t = 4)]
        threads: usize,
        #[arg(long, default_value_t = 8.0)]
        gap_open: f64,
        #[arg(long, default_value_t = 1.0)]
        gap_extend: f64,
        #[arg(long, default_value_t = 2.5)]
        edit_penalty: f64,
        /// Cost of a U-to-C ("reverse") RNA edit other than a stop
        /// read-through. `inf`, the default: none, as in seed plants,
        /// mosses, liverworts, horsetails and most lycophytes. Ferns,
        /// quillworts and hornworts edit both ways
        #[arg(long, default_value_t = f64::INFINITY)]
        reverse_edit_penalty: f64,
        /// Cost of reading a genomic in-frame TAA or TAG as Gln, or TGA as
        /// Arg, by a U-to-C edit. `inf`, the default: a stop is a stop
        #[arg(long, default_value_t = f64::INFINITY)]
        stop_edit_penalty: f64,
        /// For a profile that says where its reference proteins are edited
        /// (a 21st `edit` column): the cost of a C-to-U edit at a column
        /// edited in them. Default: --edit-penalty
        #[arg(long)]
        known_edit_penalty: Option<f64>,
        /// ...and at a column that is not. Default: --edit-penalty
        #[arg(long)]
        novel_edit_penalty: Option<f64>,
        #[arg(long, default_value_t = 300)]
        flank: i64,
        /// Cost of a frame change (skipping 1 or 2 bases) inside a gene;
        /// "inf" aligns each gene in a single reading frame
        #[arg(long, default_value_t = 20.0)]
        frameshift_penalty: f64,
        /// Profile columns each side of a frame change that must align
        /// soundly for it to be kept (default 0: keep all)
        #[arg(long, default_value_t = orfedit::FRAMESHIFT_FLANK)]
        frameshift_flank: usize,
        /// ...of which at most this many may be unaligned
        #[arg(long, default_value_t = orfedit::FRAMESHIFT_MAX_DELETED)]
        frameshift_max_deleted: usize,
        /// ...and the aligned ones must score at least this on average
        #[arg(long, default_value_t = orfedit::FRAMESHIFT_MIN_FLANK_SCORE)]
        frameshift_min_flank_score: f64,
        /// A gene keeps its frame changes only if they improve its
        /// alignment score by at least this much over a single frame
        #[arg(long, default_value_t = orfedit::FRAMESHIFT_MIN_GAIN)]
        frameshift_min_gain: f64,
        /// When no stop codon directly follows the aligned protein, look
        /// this many codons further on for one and extend the gene to it
        /// (0: only a stop right after the last aligned codon counts)
        #[arg(long, default_value_t = 50)]
        stop_search: usize,
        /// When the first aligned codon is not a start codon (ATG, or ACG
        /// an edit makes AUG), start the gene at the nearest ATG within
        /// this many codons instead (upstream, in an open frame; or
        /// downstream, at a column whose references can start with Met),
        /// else the nearest ACG (0: start where the alignment does)
        #[arg(long, default_value_t = 10)]
        start_search: usize,
        /// Cost of a stop codon inside a gene whose protein goes on
        /// matching beyond it (a pseudogene, a stop a lineage edits away,
        /// an error): the stop is aligned to a profile column, `*` in the
        /// protein, and counted in n_internal_stops. "inf": a gene ends at
        /// its first stop, or gets round it by frame changes
        #[arg(long, default_value_t = 20.0)]
        internal_stop_penalty: f64,
        /// A frame change or an internal stop with fewer than this many
        /// profile columns before or after it is not kept: the gene starts
        /// after it, or ends before it, in its own frame at its own stop
        /// (0: kept anywhere)
        #[arg(long, default_value_t = 25)]
        end_margin: usize,
        /// Genes (comma-separated; the manifest's gene name before any
        /// '|', or the profile's file name without .pssm) that start at the
        /// most upstream ATG of their open reading frame, unless they
        /// start at an ACG: for an N-terminal extension too variable to
        /// align, annotated from the frame's first ATG (angiosperm mito
        /// atp6). Looks up to 300 codons upstream
        #[arg(long, value_delimiter = ',')]
        upstream_start: Vec<String>,
    },
}

struct ManifestRow {
    species: String,
    organelle: String,
    gene: String,
    fasta: String,
    profile: String,
    contig: String,
    start: i64,
    end: i64,
    strand: char,
}

fn parse_manifest(path: &PathBuf) -> Result<Vec<ManifestRow>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading manifest {}", path.display()))?;
    let mut lines = text.lines();
    let header = lines.next().context("empty manifest")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let expected = ["species", "organelle", "gene", "fasta", "profile", "contig", "start", "end", "strand"];
    if cols != expected {
        anyhow::bail!("manifest header must be {:?}, got {:?}", expected, cols);
    }
    let mut rows = Vec::new();
    for (lineno, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 9 {
            anyhow::bail!("{}:{}: expected 9 columns, got {}", path.display(), lineno + 2, f.len());
        }
        rows.push(ManifestRow {
            species: f[0].to_string(),
            organelle: f[1].to_string(),
            gene: f[2].to_string(),
            fasta: f[3].to_string(),
            profile: f[4].to_string(),
            contig: f[5].to_string(),
            start: f[6].parse().with_context(|| format!("{}:{}: bad start", path.display(), lineno + 2))?,
            end: f[7].parse().with_context(|| format!("{}:{}: bad end", path.display(), lineno + 2))?,
            strand: f[8].chars().next().context("empty strand")?,
        });
    }
    Ok(rows)
}

struct OutRow {
    species: String,
    organelle: String,
    gene: String,
    frame: usize,
    score: f64,
    raw_start: i64,
    raw_end: i64,
    corrected_start: usize,
    corrected_end: usize,
    changed: bool,
    n_edits: usize,
    n_columns_aligned: usize,
    n_columns_total: usize,
    includes_stop: bool,
    n_internal_stops: usize,
    protein: String,
    edits: Vec<(usize, usize, char, bool)>, // genomic_pos, profile_column, resulting_aa, reverse
    frameshifts: Vec<(usize, usize, usize, usize, usize, usize, bool)>, // start, end, n_bases, profile_column, window, near_end
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::ScanBatch {
            manifest,
            gene_calls_out,
            edits_out,
            frameshifts_out,
            threads,
            gap_open,
            gap_extend,
            edit_penalty,
            reverse_edit_penalty,
            stop_edit_penalty,
            known_edit_penalty,
            novel_edit_penalty,
            flank,
            frameshift_penalty,
            frameshift_flank,
            frameshift_max_deleted,
            frameshift_min_flank_score,
            frameshift_min_gain,
            stop_search,
            start_search,
            internal_stop_penalty,
            end_margin,
            upstream_start,
        } => {
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global().ok();

            let rows = parse_manifest(&manifest)?;
            eprintln!("[info] orfedit: {} manifest rows", rows.len());

            // Preload every distinct profile once (small: ~10-20 genes total).
            let mut profile_paths: Vec<&str> = rows.iter().map(|r| r.profile.as_str()).collect();
            profile_paths.sort_unstable();
            profile_paths.dedup();
            let mut profiles: HashMap<String, Profile> = HashMap::new();
            for p in profile_paths {
                let gene = rows.iter().find(|r| r.profile == p).map(|r| r.gene.clone()).unwrap_or_default();
                match Profile::load(&PathBuf::from(p), &gene) {
                    Ok(prof) => {
                        profiles.insert(p.to_string(), prof);
                    }
                    Err(e) => eprintln!("[warn] orfedit: failed to load profile {}: {:#}", p, e),
                }
            }

            // Group rows by fasta path so each species/organelle fasta is
            // loaded exactly once, regardless of how many genes it's used for.
            let mut groups: HashMap<String, Vec<&ManifestRow>> = HashMap::new();
            for r in &rows {
                groups.entry(r.fasta.clone()).or_default().push(r);
            }
            let params = ScanParams {
                gap_open,
                gap_extend,
                edit_penalty,
                reverse_edit_penalty,
                stop_edit_penalty,
                known_edit_penalty,
                novel_edit_penalty,
                flank,
                frameshift_penalty,
                frameshift_flank,
                frameshift_max_deleted,
                frameshift_min_flank_score,
                frameshift_min_gain,
                stop_search,
                start_search,
                internal_stop_penalty,
                end_margin,
                upstream_start: false,
            };
            // the genes started at the most upstream ATG of their frame
            let upstream_genes: Vec<String> = upstream_start.iter().map(|g| g.to_lowercase()).collect();
            let upstream_for = |r: &ManifestRow| {
                let name = r.gene.split('|').next().unwrap_or("").to_lowercase();
                let stem = std::path::Path::new(&r.profile).file_stem().map(|s| s.to_string_lossy().to_lowercase()).unwrap_or_default();
                upstream_genes.iter().any(|g| *g == name || *g == stem)
            };

            let group_vec: Vec<(&String, &Vec<&ManifestRow>)> = groups.iter().collect();
            let all_results: Vec<OutRow> = group_vec
                .par_iter()
                .flat_map(|(fasta_path, group_rows)| {
                    let idx = match FastaIndex::load(&PathBuf::from(fasta_path.as_str())) {
                        Ok(idx) => idx,
                        Err(e) => {
                            eprintln!("[warn] orfedit: failed to load fasta {}: {:#}", fasta_path, e);
                            return Vec::new();
                        }
                    };
                    let mut out = Vec::with_capacity(group_rows.len());
                    for r in group_rows.iter() {
                        let profile = match profiles.get(&r.profile) {
                            Some(p) => p,
                            None => continue, // already warned above
                        };
                        let row_params = ScanParams { upstream_start: upstream_for(r), ..params };
                        match scan_gene(&idx, &r.contig, r.start, r.end, r.strand, profile, &row_params) {
                            Ok(outcome) => {
                                // "changed" = editing was actually invoked, not "coordinates differ at
                                // all" - a profile-alignment boundary almost never lands on the exact
                                // same bp as an independently-called boundary even when nothing
                                // biological changed (confirmed empirically: raw boundary noise alone
                                // made ~97% of zero-edit calls look "changed" under a coordinate-equality
                                // definition - see analysis/editing/README.md). n_edits is the actual
                                // signal this tool exists to produce.
                                let changed = !outcome.edits.is_empty();
                                out.push(OutRow {
                                    species: r.species.clone(),
                                    organelle: r.organelle.clone(),
                                    gene: r.gene.clone(),
                                    frame: outcome.frame,
                                    score: outcome.score,
                                    raw_start: r.start,
                                    raw_end: r.end,
                                    corrected_start: outcome.corrected_start,
                                    corrected_end: outcome.corrected_end,
                                    changed,
                                    n_edits: outcome.edits.len(),
                                    n_columns_aligned: outcome.n_columns_aligned,
                                    n_columns_total: outcome.n_columns_total,
                                    includes_stop: outcome.includes_stop,
                                    n_internal_stops: outcome.n_internal_stops,
                                    protein: outcome.protein,
                                    edits: outcome.edits.into_iter().map(|e| (e.genomic_pos, e.profile_column, e.resulting_aa, e.reverse)).collect(),
                                    frameshifts: outcome
                                        .frameshifts
                                        .into_iter()
                                        .map(|f| (f.start, f.end, f.n_bases, f.profile_column, f.window_start, f.window_end, f.near_end))
                                        .collect(),
                                });
                            }
                            Err(e) => {
                                eprintln!("[warn] orfedit: {} {} {}: {:#}", r.species, r.organelle, r.gene, e);
                            }
                        }
                    }
                    out
                })
                .collect();

            write_outputs(&gene_calls_out, &edits_out, frameshifts_out.as_ref(), &all_results)?;
            eprintln!(
                "[info] orfedit: {} genes scanned, {} corrected, {} with frameshifts -> {}",
                all_results.len(),
                all_results.iter().filter(|r| r.changed).count(),
                all_results.iter().filter(|r| r.frameshifts.iter().any(|f| !f.6)).count(),
                gene_calls_out.display()
            );
            Ok(())
        }
    }
}

fn write_outputs(
    gene_calls_out: &PathBuf,
    edits_out: &PathBuf,
    frameshifts_out: Option<&PathBuf>,
    results: &[OutRow],
) -> Result<()> {
    if let Some(parent) = gene_calls_out.parent() {
        fs::create_dir_all(parent)?;
    }
    if let Some(parent) = edits_out.parent() {
        fs::create_dir_all(parent)?;
    }

    let tmp_gc = gene_calls_out.with_extension("tsv.tmp");
    {
        let mut w = String::new();
        w.push_str("species\torganelle\tgene\tframe\tscore\traw_start\traw_end\tcorrected_start\tcorrected_end\tchanged\tn_edits\tn_columns_aligned\tn_columns_total\tincludes_stop\tn_frameshifts\tprotein\tn_internal_stops\n");
        for r in results {
            w.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                r.species, r.organelle, r.gene, r.frame, r.score, r.raw_start, r.raw_end,
                r.corrected_start, r.corrected_end, r.changed, r.n_edits, r.n_columns_aligned, r.n_columns_total,
                r.includes_stop, r.frameshifts.iter().filter(|f| !f.6).count(), r.protein, r.n_internal_stops
            ));
        }
        fs::write(&tmp_gc, w)?;
    }
    fs::rename(&tmp_gc, gene_calls_out)?;

    let tmp_ed = edits_out.with_extension("tsv.tmp");
    {
        let mut w = String::new();
        w.push_str("species\torganelle\tgene\tgenomic_pos\tprofile_column\tresulting_aa\tedit\n");
        for r in results {
            for (pos, col, aa, reverse) in &r.edits {
                let kind = if *reverse { "U>C" } else { "C>U" };
                w.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", r.species, r.organelle, r.gene, pos, col, aa, kind));
            }
        }
        fs::write(&tmp_ed, w)?;
    }
    fs::rename(&tmp_ed, edits_out)?;

    if let Some(path) = frameshifts_out {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let tmp = path.with_extension("tsv.tmp");
        let mut w = String::new();
        w.push_str("species\torganelle\tgene\tstart\tend\tn_bases\tprofile_column\twindow_start\twindow_end\tnear_end\n");
        for r in results {
            for (start, end, n_bases, col, ws, we, near_end) in &r.frameshifts {
                w.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    r.species, r.organelle, r.gene, start, end, n_bases, col, ws, we, near_end
                ));
            }
        }
        fs::write(&tmp, w)?;
        fs::rename(&tmp, path)?;
    }

    Ok(())
}
