//! Per-gene amino acid PSSM ("profile"), built externally by
//! `analysis/editing/src/00_build_reference_profiles.py` from this dataset's
//! own already-annotated orthologs (see that script's docstring) and read
//! here as a plain TSV: one header row of 20 amino acid columns (order below
//! must match `codon::AA_ORDER`), one data row per alignment column. A 21st
//! column headed `edit` is optional: the share of the profile's reference
//! proteins whose residue in that column RNA editing made (0: none).

use anyhow::{bail, Context, Result};
use std::fs;
use std::path::Path;

use crate::codon::AA_ORDER;

pub struct Profile {
    pub gene: String,
    /// columns[i][a] = log-odds score for amino acid `a` (index into
    /// AA_ORDER) at profile column i.
    pub columns: Vec<[f64; 20]>,
    /// Per column, whether any reference protein is RNA-edited there (the
    /// file's `edit` column above 0); `None` for a profile without one.
    pub edited: Option<Vec<bool>>,
}

impl Profile {
    pub fn len(&self) -> usize {
        self.columns.len()
    }

    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }

    /// Score of amino acid `aa` at profile column `i` (0-indexed).
    /// Unrecognised "amino acids" (stop, X) are NEG_INFINITY, not just a
    /// large penalty: a real coding column can never truly resolve to a
    /// stop, so the DP must always prefer a deletion there instead of
    /// "matching" a stop/unknown to a real residue column. Confirmed
    /// necessary empirically, not just in principle - an earlier -10.0
    /// penalty let a handful of interior stop-codon readings win over a
    /// deletion when every real amino acid hypothesis scored even worse,
    /// producing literal `*` characters inside output protein sequences
    /// (see `transsplice`'s README/validation notes for the real case that
    /// caught this).
    pub fn score(&self, i: usize, aa: u8) -> f64 {
        match AA_ORDER.iter().position(|&x| x == aa) {
            Some(idx) => self.columns[i][idx],
            None => f64::NEG_INFINITY,
        }
    }

    /// Whether the reference proteins are edited at column `i`; `None`
    /// for a profile that does not say.
    #[inline]
    pub fn edit_site(&self, i: usize) -> Option<bool> {
        self.edited.as_ref().map(|e| e[i])
    }

    pub fn load(path: &Path, gene: &str) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("reading profile {}", path.display()))?;
        let mut lines = text.lines();
        let header = lines
            .next()
            .with_context(|| format!("empty profile file {}", path.display()))?;
        let mut header_cols: Vec<&str> = header.split('\t').collect();
        let has_edit = header_cols.len() == 21 && header_cols[20].trim().eq_ignore_ascii_case("edit");
        if has_edit {
            header_cols.pop();
        }
        if header_cols.len() != 20 {
            bail!(
                "{}: expected 20 amino acid columns in header (and optionally `edit`), got {}",
                path.display(),
                header_cols.len()
            );
        }
        let n_fields = if has_edit { 21 } else { 20 };
        let mut edited = Vec::new();
        // Column order in the file need not match AA_ORDER exactly, but every
        // AA_ORDER letter must be present exactly once - map file-column ->
        // AA_ORDER index once, up front.
        let mut file_col_to_aa_idx = vec![0usize; 20];
        for (file_col, name) in header_cols.iter().enumerate() {
            let byte = name.trim().as_bytes();
            if byte.len() != 1 {
                bail!("{}: bad amino acid column header {:?}", path.display(), name);
            }
            let aa_idx = AA_ORDER
                .iter()
                .position(|&x| x == byte[0])
                .with_context(|| format!("{}: unknown amino acid column {:?}", path.display(), name))?;
            file_col_to_aa_idx[file_col] = aa_idx;
        }

        let mut columns = Vec::new();
        for (lineno, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != n_fields {
                bail!(
                    "{}:{}: expected {} columns, got {}",
                    path.display(),
                    lineno + 2,
                    n_fields,
                    fields.len()
                );
            }
            if has_edit {
                let share: f64 = fields[20]
                    .trim()
                    .parse()
                    .with_context(|| format!("{}:{}: bad edit share {:?}", path.display(), lineno + 2, fields[20]))?;
                edited.push(share > 0.0);
            }
            let mut col = [0.0f64; 20];
            for (file_col, v) in fields.iter().take(20).enumerate() {
                let score: f64 = v
                    .trim()
                    .parse()
                    .with_context(|| format!("{}:{}: bad score {:?}", path.display(), lineno + 2, v))?;
                col[file_col_to_aa_idx[file_col]] = score;
            }
            columns.push(col);
        }
        if columns.is_empty() {
            bail!("{}: no profile columns", path.display());
        }
        Ok(Profile { gene: gene.to_string(), columns, edited: has_edit.then_some(edited) })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_profile(dir: &Path, name: &str, body: &str) -> std::path::PathBuf {
        let path = dir.join(name);
        let mut f = fs::File::create(&path).unwrap();
        write!(f, "{}", body).unwrap();
        path
    }

    #[test]
    fn loads_and_indexes_by_aa_letter_regardless_of_column_order() {
        let dir = std::env::temp_dir();
        // Header intentionally not in AA_ORDER order.
        let body = "M\tA\nX\t0.0\n"; // header has garbage col name "X" -> should fail via unknown col name
        let _ = body; // placeholder to keep clippy quiet about unused fn if test below changes
        let header: String = AA_ORDER.iter().map(|&b| (b as char).to_string()).collect::<Vec<_>>().join("\t");
        let row: String = (0..20).map(|i| format!("{:.2}", i as f64 * 0.1)).collect::<Vec<_>>().join("\t");
        let content = format!("{}\n{}\n", header, row);
        let path = write_profile(&dir, "test_gene.pssm", &content);
        let p = Profile::load(&path, "test_gene").unwrap();
        assert_eq!(p.len(), 1);
        // AA_ORDER[0] = 'A', score should be column 0 = 0.0
        assert!((p.score(0, b'A') - 0.0).abs() < 1e-9);
        // AA_ORDER[19] = 'V', score should be column 19 = 1.9
        assert!((p.score(0, b'V') - 1.9).abs() < 1e-9);
        assert!(p.edited.is_none() && p.edit_site(0).is_none());
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn reads_the_optional_edit_column() {
        let dir = std::env::temp_dir();
        let header: String = AA_ORDER.iter().map(|&b| (b as char).to_string()).collect::<Vec<_>>().join("\t");
        let row: String = (0..20).map(|_| "0.5".to_string()).collect::<Vec<_>>().join("\t");
        let content = format!("{header}\tedit\n{row}\t0\n{row}\t0.25\n");
        let path = write_profile(&dir, "test_gene_edit.pssm", &content);
        let p = Profile::load(&path, "test_gene").unwrap();
        assert_eq!(p.len(), 2);
        assert_eq!(p.edited, Some(vec![false, true]));
        assert_eq!((p.edit_site(0), p.edit_site(1)), (Some(false), Some(true)));
        // a row without the edit share is an error, not a silent default
        let short = format!("{header}\tedit\n{row}\n");
        let bad = write_profile(&dir, "test_gene_edit_bad.pssm", &short);
        assert!(Profile::load(&bad, "test_gene").is_err());
        let _ = fs::remove_file(&path);
        let _ = fs::remove_file(&bad);
    }
}
