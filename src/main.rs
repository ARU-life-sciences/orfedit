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
        #[arg(long, default_value_t = 4)]
        threads: usize,
        #[arg(long, default_value_t = 8.0)]
        gap_open: f64,
        #[arg(long, default_value_t = 1.0)]
        gap_extend: f64,
        #[arg(long, default_value_t = 2.5)]
        edit_penalty: f64,
        #[arg(long, default_value_t = 300)]
        flank: i64,
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
    protein: String,
    edits: Vec<(usize, usize, char)>, // genomic_pos, profile_column, resulting_aa
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Command::ScanBatch { manifest, gene_calls_out, edits_out, threads, gap_open, gap_extend, edit_penalty, flank } => {
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
            let params = ScanParams { gap_open, gap_extend, edit_penalty, flank };

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
                        match scan_gene(&idx, &r.contig, r.start, r.end, r.strand, profile, &params) {
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
                                    protein: outcome.protein,
                                    edits: outcome.edits.into_iter().map(|e| (e.genomic_pos, e.profile_column, e.resulting_aa)).collect(),
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

            write_outputs(&gene_calls_out, &edits_out, &all_results)?;
            eprintln!(
                "[info] orfedit: {} genes scanned, {} corrected -> {}",
                all_results.len(),
                all_results.iter().filter(|r| r.changed).count(),
                gene_calls_out.display()
            );
            Ok(())
        }
    }
}

fn write_outputs(gene_calls_out: &PathBuf, edits_out: &PathBuf, results: &[OutRow]) -> Result<()> {
    if let Some(parent) = gene_calls_out.parent() {
        fs::create_dir_all(parent)?;
    }
    if let Some(parent) = edits_out.parent() {
        fs::create_dir_all(parent)?;
    }

    let tmp_gc = gene_calls_out.with_extension("tsv.tmp");
    {
        let mut w = String::new();
        w.push_str("species\torganelle\tgene\tframe\tscore\traw_start\traw_end\tcorrected_start\tcorrected_end\tchanged\tn_edits\tn_columns_aligned\tn_columns_total\tincludes_stop\tprotein\n");
        for r in results {
            w.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                r.species, r.organelle, r.gene, r.frame, r.score, r.raw_start, r.raw_end,
                r.corrected_start, r.corrected_end, r.changed, r.n_edits, r.n_columns_aligned, r.n_columns_total,
                r.includes_stop, r.protein
            ));
        }
        fs::write(&tmp_gc, w)?;
    }
    fs::rename(&tmp_gc, gene_calls_out)?;

    let tmp_ed = edits_out.with_extension("tsv.tmp");
    {
        let mut w = String::new();
        w.push_str("species\torganelle\tgene\tgenomic_pos\tprofile_column\tresulting_aa\n");
        for r in results {
            for (pos, col, aa) in &r.edits {
                w.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\n", r.species, r.organelle, r.gene, pos, col, aa));
            }
        }
        fs::write(&tmp_ed, w)?;
    }
    fs::rename(&tmp_ed, edits_out)?;

    Ok(())
}
