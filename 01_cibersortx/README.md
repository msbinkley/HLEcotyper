# 1. CIBERSORTx — signature matrix, benchmarking and deconvolution

All jobs run CIBERSORTx from the official Singularity images (`fractions.sif`,
`hires.sif`; obtain from https://cibersortx.stanford.edu after registering) on a SLURM
cluster. Every script binds one input directory to `/src/data` inside the container, so
**all `/src/data/...` paths below refer to files placed in `$INPUT_DIR`** (default
`$CIBERSORT_ROOT/bulk_input`). Set once per shell:

```bash
export CIBERSORTX_USER="you@institution.edu"     # registered CIBERSORTx account
export CIBERSORTX_TOKEN="<token from the website>"
export CIBERSORT_ROOT=/path/to/Cibersort         # holds *.sif and bulk_input/
```

| step | script | mode | input (in `$INPUT_DIR`) | key output |
|---|---|---|---|---|
| build | `01_build_signature_matrix.sh` | Fractions, `--single_cell TRUE --fraction 0 --filter TRUE --G.min 400 --G.max 800` | `cHL_scRNA_reference_4000.txt` (genes × 44,769 cells, header = cell type) | signature matrix (16 cell types) → `resources/cHL_signature_matrix_400_800.txt` |
| benchmark 1 | `02_benchmark_pseudobulk.sh` | Fractions | `pseudobulk_cHL_<celltype>_cpm.txt` per cell type (also for external 10X / Aoki / Stewart references) | recovered vs known fraction |
| benchmark 1b | `03_benchmark_HRS_subsample.sh` | Fractions (rebuild) | `SubSample_HRS_<p>pct_rep1.txt` | signature-matrix stability vs HRS cell count |
| **deconvolution** | `04_hires_bulk_tumor.sh` | HiRes, `--rmbatchSmode TRUE --subsetgenes coding_gene_list.txt` | `cHL_bulk_TPM.txt` (187 tumors) | `CIBERSORTxGEP_NA_Fractions-Adjusted.txt` (Fig. 1b) + `CIBERSORTxHiRes_NA_<celltype>_Window*.txt` (EcoTyper input) |
| Visium | `05_hires_visium.sh` | HiRes S-mode | `ST_data_<tumor>.txt` ×4 | per-spot fractions → EcoTyper Visium recovery |
| cfRNA | `06_hires_cfRNA.sh` | HiRes S-mode | `cfRNA_coding_symbol_TPM.txt` | cfRNA HRS fraction; EcoTyper cfRNA recovery |

Benchmarks 2 (VAF) and 3 (PhenoCycler) compare the HRS fraction from step 4 against
mutation-based purity and imaging counts; they are plotting/statistics only and need no
CIBERSORTx run.

## Inputs you must provide

* **Bulk TPM** (`cHL_bulk_TPM.txt`): tab-delimited, first column `Gene` (HGNC symbol),
  one column per sample, restricted to protein-coding genes (`resources/coding_gene_list.txt`)
  and renormalised so each column sums to 1e6. How the TPM is produced upstream is out of scope
  here — any Salmon/RSEM gene-level TPM works.
* **Single-cell reference** (`cHL_scRNA_reference_4000.txt`): counts for the 4,000 variable
  features of the annotated scRNA/snRNA-seq object; header row gives each cell's label
  (16 labels listed in `resources/cell_type_mapping_16_to_13.tsv`). Exported with Seurat
  `GetAssayData(obj, slot = "counts")[VariableFeatures(obj), ]` and `Idents(obj)` as column names.

## Resources shipped

* `resources/cHL_signature_matrix_400_800.txt` — the custom cHL signature matrix (1,845 genes × 16 cell types).
  Column 1 is `NAME`; gene symbols are in the namespace of the reference used at build time
  (GENCODE v27-era symbols; CIBERSORTx applies `make.names()` internally so hyphens appear as dots).
* `resources/coding_gene_list.txt` — protein-coding genes passed to `--subsetgenes`.
* `resources/cell_type_mapping_16_to_13.tsv` and `resources/pool_fractions_16_to_13.R` —
  the a-priori pooling of 16 → 13 cell types (`B_preGC + B_Memory → B_cell`, `NK + ILC → NK_cell`,
  `T_MKI67` dropped, rows renormalised) that produces the fraction file EcoTyper consumes.

## Notes

* Runtime: signature-matrix build ≈ hours (10 CPU); bulk HiRes ≈ 1–2 h on 14 CPUs.
* CIBERSORTx HiRes writes `NA` for genes it cannot impute and a constant value for genes it
  cannot resolve per sample; treat both as missing downstream.
