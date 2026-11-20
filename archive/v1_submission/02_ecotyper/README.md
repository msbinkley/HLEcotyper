# 2. EcoTyper — cell-state and Hodgkin lymphoma ecotype (HLE) discovery and recovery

Uses the unmodified EcoTyper framework (https://github.com/digitalcytometry/ecotyper,
R 4.2). Clone it, then run the scripts here **from inside that clone** — EcoTyper resolves
`-d <Discovery dataset name>` against its own `EcoTyper/<name>/` folder.

```bash
git clone https://github.com/digitalcytometry/ecotyper ecotyper-master
cd ecotyper-master
cp /path/to/HLEcotyper/02_ecotyper/discovery_config.yaml .
sbatch /path/to/HLEcotyper/02_ecotyper/01_discovery.sbatch
```

## Discovery (`discovery_config.yaml`, `01_discovery.sbatch`)

| parameter | value | Methods wording |
|---|---|---|
| Expression matrix | bulk TPM, 187 tumors | same file as `01_cibersortx/04_hires_bulk_tumor.sh` |
| Cell type fractions | 13 pooled fractions (`pool_fractions_16_to_13.R`) | "pooling the 16 cell types into 13" |
| cell-type-specific expression | CIBERSORTx HiRes, S-mode, coding genes | "CIBERSORTx HiRes module" |
| NMF genes per cell type | top 1,000 most variable (EcoTyper default) | |
| Number of NMF restarts | 50 | |
| Maximum number of states per cell type | 6 | |
| Cophenetic coefficient cutoff | 0.90 | |
| Minimum number of states in ecotypes | 6 | |
| Pipeline steps to skip | `[7]` (full run) — published model: `[1,2,3,7]` on the pre-computed HiRes output | |
| ecotype clustering | hypergeometric P < 0.01 co-occurrence → Jaccard → hierarchical, ≥ 4 states per cluster | EcoTyper defaults |

`[1,2,3,7]` skips fraction estimation, HiRes purification and NMF, and re-runs rank
selection → ecotype discovery on an existing `EcoTyper/<name>/` tree; that is how the
cophenetic sweep (`07_cophenetic_sweep.sh`) and the final model were produced after the
HiRes run of `01_cibersortx/04_hires_bulk_tumor.sh`.

### Random seeds
EcoTyper fixes its seeds internally (nothing to set in the config):
`pipeline/state_discovery_NMF.R` uses `seed = 1234 + restart` for each of the 50 NMF
restarts (`method = "brunet"`); `state_discovery_*_filter_genes.R` and
`state_recovery_scRNA.R` use `set.seed(1234)`; `lib/ecotyper.R` uses `set.seed(100 + x)`
for the ecotype clustering bootstrap and `lib/misc.R` `set.seed(0)`.

### Trained model (`model/`)
`model/` holds the text outputs of the final discovery run
(`Output folder = cHL_output_90_6_6`, 2 ecotypes):

* `rank_data.txt` — cophenetic curve / chosen number of states per cell type
* `<cell type>/gene_info.txt` — **cell-state-defining genes** (state, max fold-change, initial state) for the 13 cell types
* `<cell type>/state_assignment.txt`, `state_abundances.txt`, `heatmap_*.txt` — per-sample states
* `Ecotypes/ecotypes.txt` — state → HLE membership (E1 = HLE1, E2 = HLE2); `ecotype_assignment.txt`, `ecotype_abundance.txt` — per-sample HLE calls; `discovery/` — Jaccard clustering intermediates (`mapping_all_states.txt`, `assignment_p_vals.txt`, `silhouette*.txt`)

The *recovery-ready* model is the `EcoTyper/cHL_discovery/` tree that EcoTyper writes inside
its clone (NMF W matrices, purified profiles; ~GB-scale). It is too large for GitHub and is
deposited on Zenodo with the paper; unpack it into your EcoTyper clone as
`ecotyper-master/EcoTyper/cHL_discovery/` before running any recovery script below.

## Recovery

| script | data | EcoTyper mode | notes |
|---|---|---|---|
| `02_recovery_bulk_microarray.sh` | Steidl 2010 microarray, 130 cases (GSE17920) | `EcoTyper_recovery_bulk.R -a annotation -c Age45,Age60` | matrix batch-corrected to the signature matrix with CIBERSORTx B-mode first |
| `03_recovery_bulk_cfRNA.sh` | plasma cfRNA (TMMwsp-CPM, log2, RUVg) | `EcoTyper_recovery_bulk.R` | |
| `04_recovery_bulk_EPICseq.sh` | EPIC-seq inferred expression (linear) | `EcoTyper_recovery_bulk.R` | |
| `05_recovery_scRNA_external.sh` | 3 external cHL scRNA-seq sets (34 cases) + in-house scRNA of the Visium tumors | `EcoTyper_recovery_scRNA.R -c Dataset -s 10000` | labels must use the 13 pooled cell-type names |
| `06_recovery_visium.sh` + `visium_recovery_config_template.yml` | 4 Visium tumors | `EcoTyper_recovery_visium.R` | per-spot fractions from CIBERSORTx (`*_Ciber`) or CytoSPACE (`*_Cyto`) |
| `07_cophenetic_sweep.sh` | — | steps 4–8 at cutoffs 0.90–1.00 + microarray recovery | sensitivity analysis |

Outputs of every recovery run mirror the discovery layout
(`<output>/<dataset>/<cell type>/state_abundances.txt`, `<output>/<dataset>/Ecotypes/ecotype_abundance.txt`).
Differences in state abundance between HLEs were tested with two-sided Welch t-tests and
Q-values (Supplementary Table S18) on `state_abundances.txt` — plain R, no script needed.
