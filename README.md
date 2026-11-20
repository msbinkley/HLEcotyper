# HLEcotyper

Code and resources for **"Distinct cell state ecosystems for classic Hodgkin lymphoma"**
(Su, Subramanian et al.; M.S. Binkley, corresponding author).

The study builds a cHL-specific single-cell reference (110,043 cells, 16 cell types) and CIBERSORTx signature
matrix, deconvolutes bulk RNA-seq of 187 tumors, discovers 30 cell states across 13 cell types and two Hodgkin
lymphoma ecotypes (HLE1, HLE2) with EcoTyper, and recovers them in an external microarray cohort (GSE17920),
the single-cell data, Visium spatial transcriptomics, plasma cfRNA, EPIC-seq and an ex vivo tumor culture.

This repository holds the analysis scripts **as they were run** for the published model
(EcoTyper model `r16_c0.933_m4_e6`) and the small resources needed to reuse it. Scripts that only draw figures
are not included. Processed data are in the paper's Supplementary Data (expression tables and Seurat objects).

## Layout

| folder | what it does | manuscript |
|---|---|---|
| [`01_scRNA_reference/`](01_scRNA_reference/) | merge of the five single-cell datasets on one gene universe, Harmony integration, cluster annotation from published labels, lineage sub-UMAP clean-up, unstable-cell loops; ends in the round-16 reference (`v2_C22_round16_recluster_and_check.R`); marker genes (`v3_J1`) | Methods "scRNA-seq annotation"; Fig. 2a, Fig. S1 |
| [`02_cibersortx/`](02_cibersortx/) | Ensembl/symbol clean-up of the signature matrix (`v2_D_undot_sigmatrix.py`), 16 → 13 pooling (`v2_D_pool13.R`), microarray symbols (`v2_E1`), plasma cfRNA fractions (`v3_I6a`, `v3_I6b`), LOD95 spike-in (`v3_I10a`, `v3_I13a`, `v3_I13b`) | Methods "Construction and benchmarking…", "Cell-free RNA…", "LOD95…"; Fig. 2, 6b–d |
| [`03_ecotyper/`](03_ecotyper/) | round-16 chain: reference → CIBERSORTx signature matrix and S-mode fractions → EcoTyper discovery → model search with microarray recovery (`v2_H1`, `v2_H2`, `v2_D11`, `v2_D14`); the final model and its recoveries (`v2_H4_model_chain.sh r16_c0.933_m4_e6`, `v2_H4_bulk_finish.sh`, `v2_H4_scrna_recovery_inputs.R`); Visium recovery (`v3_I14`); ex vivo recovery (`v2_G5`, `v2_G6`) | Methods "Discovery of cell states…", "Validation…"; Fig. 3–7 |
| [`04_cytospace/`](04_cytospace/) | CytoSPACE inputs from the round-16 reference and runs on the four Visium tumors | Methods "Spatial transcriptomics"; Fig. 5 |
| [`05_trust4/`](05_trust4/) | TRUST4 (v1.1.10) on fastp-trimmed bulk RNA-seq, TCR/BCR Shannon entropy | Methods "Bulk tumor RNA sequencing"; Fig. S4a |
| [`06_microbiome/`](06_microbiome/) | Kraken2/Bracken pipeline and EBV / microbiome statistics (unchanged from the first release) | Methods "Microbial pathogen analysis"; Fig. S7 |
| [`resources/`](resources/) | round-16 signature matrix, 16 → 13 mapping, final EcoTyper model | — |
| [`archive/v1_submission/`](archive/v1_submission/) | the first release (code and model of the originally submitted version) | — |

Every script starts with a header that says what it does, its inputs and outputs, how it was run and which random
seeds it uses.

## Resources

* `resources/cHL_signature_matrix_round16.txt` — the cHL CIBERSORTx signature matrix (1,698 genes × 16 cell types)
* `resources/cell_type_mapping_16_to_13.tsv` — pooling of the 16 deconvoluted cell types into the 13 EcoTyper cell
  types (B_GC + B_Memory → B_postGC; NK + ILC → NK; Cycling excluded)
* `resources/ecotyper_model/` — the final model `r16_c0.933_m4_e6` (cophenetic cutoff 0.933, at most 4 states per
  cell type, at least 6 states per ecotype): per cell type the state-defining genes (`gene_info.txt`), state
  abundances and assignments of the 187 training tumors and EcoTyper's heat-map data; `Ecotypes/` holds the
  20 HLE-defining states (`ecotypes.txt`, E1 = HLE1, E2 = HLE2) and the per-tumor HLE abundance and assignment

## Running the scripts

The scripts were run on a Slurm cluster (Stanford Sherlock) with R 4.2.0, and keep the folder layout of that run.
Paths are written with these variables; set them to your own folders:

| variable | meaning |
|---|---|
| `PROJECT_ROOT` | the project folder that holds `cHL_v2/`, `Ecotyper/`, `Cibersort/`, `Cytospace/` |
| `LAB_ROOT` | shared lab storage (sequencing data, shared references) |
| `BOX_ROOT` | the local project folder used for the pulls to the laptop |
| `ANNOTATION_DIR` | gene annotation files (GENCODE / Ensembl GTFs) |
| `MAIL_USER` | Slurm e-mail address for job notices |

**CIBERSORTx** needs a registered account (https://cibersortx.stanford.edu) and its Singularity images. Scripts that
call it source `v2_D_csx_credentials.sh`, which is not part of this repository: write it yourself, setting
`CSX_USER` and `CSX_TOKEN`, and keep it out of version control.
**EcoTyper** (https://github.com/digitalcytometry/ecotyper) scripts are run from inside an EcoTyper clone.
**CytoSPACE** v1.1.0 (https://github.com/digitalcytometry/cytospace). **TRUST4** v1.1.10.

Sample identifiers in the scripts are the de-identified codes used in the Supplementary Data and the released Seurat
objects (for example `tumor1`–`tumor4` for the Visium tumors).

## Citation and license

Please cite the paper and the tools used (CIBERSORTx, EcoTyper, CytoSPACE, Harmony, Seurat, TRUST4, Kraken2,
Bracken). License: see `LICENSE`.
