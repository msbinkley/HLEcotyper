# 5. Intratumoral microbiome and EBV (Fig. 7, Supplementary Fig. S7, Table S19)

Two parts: **`pipeline/`** turns FFPE bulk RNA-seq FASTQs into per-sample species calls
(Kraken2 → Bracken → filtering), and **`analysis/`** performs the EBV-centred statistics and
figures on those calls.

## 5a. Profiling pipeline (`pipeline/`, SLURM)

Host depletion (STAR to GRCh38, rRNA/low-complexity scrub, T2T-CHM13 second pass keeping
only read pairs with both mates unmapped) → Kraken2 (PlusPF database, human genome included,
`--confidence 0.1 --minimum-hit-groups 2`, minimizer report) → Bracken species re-estimation
(read length 100, threshold 5) → depth normalisation to reads per million input read pairs (RPM)
→ filters (reads ≥ 5, RPM ≥ 0.5, prevalence ≥ 1 %, distinct minimizers ≥ 5, minimizer coverage
≥ 0.01, host and reagent/kitome blocklist removed). All knobs live in `pipeline/config.sh`;
`pipeline/README.md` documents every step and output.

```bash
cd pipeline
# edit config.sh (database, container, reference paths) — or override in the environment:
TRIM_DIR=/path/to/cHL_fastp \
RUN_DIR=/path/to/cHL_microbiome \
SAMPLE_MANIFEST= \
bash run_all.sh 8              # submits 01 -> 02 -> 03 -> 05 as a SLURM afterok chain
# final tables: $RUN_DIR/05_results/cross_validated_species.csv (all calls, RPM, flags)
#               $RUN_DIR/05_results/high_confidence_taxa.csv
```

Containers/databases are built once with `build_containers.sh`, `download_pluspf_db.sh`,
`build_human_rrna.sh`, `build_t2t_index.sh`; `preflight.sh` checks the setup.

## 5b. Downstream analysis (`analysis/`, R ≥ 4.2)

Run from `analysis/` with inputs in `analysis/data/`:

| file in `data/` | source |
|---|---|
| `cHL_cross_validated_species.csv` | `05_results/cross_validated_species.csv` from the pipeline |
| `sample_design_clean.csv` | `sample, cohort, ebv (pos/neg/NA), keep` — EBV = clinical EBER/IHC; duplicates and mis-classified samples flagged `keep = FALSE`; analyses restricted to 167 tumors (75 EBV+, 92 EBV−) |
| `ebv_concordance_input.csv` | per-sample `sample, clinical_ebv, ebv_reads, ebv_rpm, ebv_high_conf, ebv_present` for HHV-4 (taxid 3050299 / 10376) |

```bash
Rscript diff_abundance.R                 # Wilcoxon-on-RPM + Fisher, prevalence >= 10 %, BH-FDR; per-comparison CSVs
Rscript plot_volcano_broken.R            # 4-panel volcano (broken axis for EBV)
Rscript plot_volcano_EBV_cHL.R           # EBV+ vs EBV- cHL volcano, EBV excluded  (Fig. 7d)
Rscript plot_ebv_positive_control_cHL.R  # EBV RPM by EBV status, Wilcoxon          (Fig. 7a)
Rscript ebv_concordance.R                # sensitivity/specificity/kappa; AUC with DeLong CI
Rscript plot_roc_youden.R                # ROC of EBV RPM vs clinical status          (Fig. 7b)
Rscript landscape_pcoa_alpha.R           # Robust-Aitchison PCoA, PERMANOVA (adonis2, 9,999 perms),
                                         # PERMDISP (betadisper), with and without EBV (Fig. 7c)
```

Outputs go to `table/` (CSV) and `analysis/` (PDF). Contaminant genera are flagged, not
removed, in the per-taxon tables (`contaminant` column) — see `METHOD.md` for the full
rationale, caveats (low-biomass FFPE, cross-cohort batch confounding) and control taxa.
All scripts call `set.seed(1)`.
