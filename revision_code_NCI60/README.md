# NCI60 revision code

This directory is a portable, code-only extraction of the NCI60 analyses used
for the manuscript revision. Copy the complete `revision_code_NCI60` directory
to the root of the NCI60 analysis repository before committing it there.

It contains no input data, cached results, or generated outputs. By default,
each script writes into `revision_code_NCI60/outputs/`, which is ignored by the
included `.gitignore`. Set `REVISION_OUTPUT_DIR` to write elsewhere.

## Included analyses

| Script | Revision analysis | Main outputs |
| --- | --- | --- |
| `R1.4_nci60_metadata_factor_fdr.R` | BH correction for NCI60 metadata-to-MOFA-factor associations | `R1.4_nci60_metadata_factor_fdr.csv`, `R1.4_nci60_tissue_factor_fdr.csv` |
| `R1.5_nci60_threshold_sensitivity.R` | Candidate-threshold sensitivity for the Factor 4 TF-to-ligand MOON branch | candidate, score, and pathway-control stability tables |

## Run

From the NCI60 repository root:

```sh
Rscript revision_code_NCI60/R1.4_nci60_metadata_factor_fdr.R
Rscript revision_code_NCI60/R1.5_nci60_threshold_sensitivity.R
```

The scripts infer the repository root as the parent of this directory. If the
folder is kept elsewhere, set `NCI60_REPO_DIR` to the repository root.

## Requirements

R1.4 uses only base R and expects:

```text
results/mofa/z_matrix.RData       (object: Z_matrix)
support/all_metadata.RData        (object: all_metadata)
```

R1.5 expects the following inputs in their existing repository locations:

```text
results/mofa/mofa_weights.RData
data/cosmos/ligrec_TF_moon_inputs.Rdata
support/dorothea_df.RData
data/RNA/RNA_log2_FPKM_clean.csv
results/cosmos/moon/moon_res_rec_to_TFmet.csv
results/cosmos/moon/meta_network_filtered.csv
support/c2.cp.v2022.1.Hs.symbols.gmt
```

It requires `cosmosR`, `readr`, `dplyr`, and `stringr`. The analysis calls a
small number of `cosmosR` internal functions; it was developed against
`cosmosR` 1.19.1, so use a compatible version for an exact re-run.

## Scope boundary

The R1.4 source was deliberately rewritten as an NCI60-only analysis. This
package contains no other cohort-analysis code or outputs.
