# HLEcotyper

Code and resources for **"Distinct cell state ecosystems for classic Hodgkin lymphoma"**
(Su\*, Subramanian\* et al.; M.S. Binkley, corresponding). The study builds a cHL-specific
CIBERSORTx signature matrix from single-cell/single-nucleus RNA-seq, deconvolutes 187 bulk
tumors, discovers 32 cell states and two Hodgkin lymphoma ecotypes (HLE1, HLE2) with
EcoTyper, recovers them in external cohorts, spatial transcriptomics, plasma cfRNA and
EPIC-seq, and profiles the intratumoral microbiome/EBV.

This repository documents **how each analysis was run** — the exact commands, parameters
and file conventions — together with the small resources needed to reuse the model. Raw and
processed data are distributed as described in the paper's Data Availability statement,
and each pipeline starts from files you already have (e.g. a gene-level TPM matrix).

## Layout

| folder | pipeline | manuscript |
|---|---|---|
| [`01_cibersortx/`](01_cibersortx/) | signature-matrix construction, benchmarking, HiRes deconvolution of bulk tumors, Visium spots and cfRNA | Fig. 1–2, 5–6; Methods "Construction and benchmarking…", "Deconvolution…" |
| [`02_ecotyper/`](02_ecotyper/) | cell-state / ecotype **discovery** and **recovery** (microarray, cfRNA, EPIC-seq, external scRNA-seq, Visium) + trained-model outputs | Fig. 3–6; Methods "Discovery of cell states…", "Validation…" |
| [`03_cytospace/`](03_cytospace/) | single-cell → Visium spot assignment | Fig. 5; Methods "Spatial transcriptomics" |
| [`04_cfRNA_LOD95/`](04_cfRNA_LOD95/) | in-silico HRS spike-in and LOD95 logistic model | Fig. 6c–d; Methods "LOD95 analysis…" |
| [`05_microbiome/`](05_microbiome/) | Kraken2/Bracken profiling pipeline + EBV statistics (Wilcoxon, ROC, PCoA/PERMANOVA, differential abundance) | Fig. 7; Methods "Microbial pathogen analysis" |

Each folder has its own README with a run table; scripts are SLURM submission wrappers
(`sbatch`) or R scripts and are annotated with the inputs they expect.

## Prerequisites

* **CIBERSORTx** Singularity images (`fractions.sif`, `hires.sif`) and a registered account —
  https://cibersortx.stanford.edu. Export `CIBERSORTX_USER` and `CIBERSORTX_TOKEN`; scripts
  refuse to run without them and no credentials are stored in this repository.
* **EcoTyper** — https://github.com/digitalcytometry/ecotyper (R 4.2); the scripts in
  `02_ecotyper/` are run from inside that clone.
* **CytoSPACE** v1.1 — https://github.com/digitalcytometry/cytospace (`03_cytospace/install.sh`).
* **Kraken2 2.1.x / Bracken 3**, PlusPF database, STAR, bowtie2, BBMap — see `05_microbiome/pipeline/README.md`.
* R packages: `Seurat`, `dplyr`, `ggplot2`, `vegan`, `pROC`, `patchwork`, `ggrepel`, `edgeR`, `RUVSeq`.
* SLURM (`--partition`, memory and CPU requests are set at the top of each script and are easy to change).

## Inputs assumed to exist

| input | format | used by |
|---|---|---|
| bulk tumor TPM | tab-delimited, `Gene` + one column per sample, protein-coding genes, columns sum to 1e6 | 01 (deconvolution), 02 (discovery) |
| single-cell reference | counts for 4,000 variable features × 44,769 annotated cells; header = cell-type label | 01 (signature matrix, S-mode reference), 03 |
| Space Ranger output per Visium tumor | standard `outs/` | 01 (Visium HiRes), 02 (Visium recovery), 03 |
| cfRNA TPM; EPIC-seq inferred expression; external scRNA-seq counts; microarray matrix (GSE17920) | genes × samples/cells | 01, 02, 04 |
| fastp-trimmed FASTQs | `R1_<sample>.trimmed.fastq.gz` | 05 |

## Resources provided

* `01_cibersortx/resources/cHL_signature_matrix_400_800.txt` — the custom cHL signature matrix (16 cell types)
* `01_cibersortx/resources/coding_gene_list.txt` — protein-coding gene list used with `--subsetgenes`
* `01_cibersortx/resources/cell_type_mapping_16_to_13.tsv` + `pool_fractions_16_to_13.R` — the 16 → 13 cell-type pooling
* `02_ecotyper/model/` — trained discovery outputs: number of states per cell type, cell-state-defining genes (`<cell type>/gene_info.txt`), per-sample state and ecotype assignments, state → HLE membership
* `02_ecotyper/README.md` — random seeds used by EcoTyper and the recovery-mode commands
* `05_microbiome/pipeline/contaminants_blocklist.txt` — reagent/kitome blocklist

The full recovery-ready EcoTyper model (`EcoTyper/cHL_discovery/`, GB-scale) is deposited
on Zenodo alongside the paper and is dropped into an EcoTyper clone before running recovery.

## Reproducing the main figures

The figure panels derive from the tables produced above: cell fractions
(`CIBERSORTxGEP_NA_Fractions-Adjusted.txt`), state abundances/assignments and ecotype
abundances (`02_ecotyper/model/…`), CytoSPACE spot assignments, the LOD95 fit
(`04_cfRNA_LOD95/03_fit_LOD95.R`) and the microbiome statistics (`05_microbiome/analysis/*.R`,
which write their figure PDFs directly). Survival analyses (Kaplan–Meier, Cox) were run in R
(`survival`, `survminer`) on the per-sample HLE assignment and clinical follow-up in the
Supplementary Data.

## Related tools

TRUST4 (immune-repertoire reconstruction from the same bulk RNA-seq; v1.0.9 container,
`run-trust4 -f hg38_bcrtcr.fa --ref human_IMGT+C.fa -1 R1 -2 R2` per sample) was used for
the repertoire analyses referenced in the paper and is not part of the five pipelines above.

## Citation and license

Please cite the paper and the underlying tools (CIBERSORTx, EcoTyper, CytoSPACE, Kraken2,
Bracken). Released under the MIT License (see `LICENSE`); a citable Zenodo DOI is minted on
publication.
