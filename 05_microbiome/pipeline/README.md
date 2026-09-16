# Tumor-microbiome pipeline (FFPE RNA-seq)

Detect and quantify microbial reads from already QC/host-removed FFPE RNA-seq
FASTQs. The default workflow starts from files named:

```text
${sample}_microbe_R1.fastq
${sample}_microbe_R2.fastq
```

Run cHL and NLPHL cohorts separately by overriding `COHORT_DIR`, `D_HOST`,
`RUN_DIR`, and/or `SAMPLE_LIST` in the environment. Step 01 host/T2T/rRNA cleanup
scripts are retained for manual reruns, but `run_all.sh` no longer submits them.

Final calls are single-classifier Kraken2/Bracken calls with KrakenUniq-style
minimizer, prevalence, abundance, host, and blocklist filtering.

---

## Workflow at a glance

`workflow.mermaid` renders the same default flow:

```text
Host-removed microbe FASTQs
        |
STEP 02 Kraken2
  confidence 0, min-hit-groups 2, minimizer report
        |
STEP 03 Bracken species
  READ_LEN 100, threshold 10, species BIOM, minimizer table, QC
        |
STEP 05 aggregate + filter
  remove host
  flag/drop blocklist contaminants
  require distinct minimizers >= 5 and minimizer coverage >= 0.01
  require reads >= 5 and RPM >= 0.5
  require prevalence >= 0.01 of cohort
        |
high_confidence_taxa.csv + cohort_high_confidence_summary.csv
        |
STEP 06 figures + pipeline_summary.txt
```


---

## How to run

```bash
cd script/taxonomic/pipeline
nano config.sh
bash run_all.sh 8
squeue -u $USER
```

`run_all.sh` reuses `$SAMPLE_LIST` if it already exists. If not, it builds one
from `$D_HOST/*_microbe_R1.fastq` and matching `_R2.fastq` files.

Final tables land in `$RUN_DIR/05_results/`; figures land in `$RUN_DIR/06_figures/`.

Manual submission:

```bash
sbatch --array=1-N%8 02_kraken2.sbatch
sbatch --dependency=afterok:<job02> 03_bracken_merge.sbatch
sbatch --dependency=afterok:<job03> 05_downstream.sbatch
```

---

## Steps

| step | script | default status | what it does |
|------|--------|----------------|--------------|
| 00 | `00_build_sample_sheet.sh` | manual/legacy | builds a sample sheet from trimmed FASTQs for step 01-era runs |
| 01 | `01_host_removal.sbatch` | already completed upstream | STAR host removal plus rRNA, low-complexity, dedup, optional T2T cleanup |
| 02 | `02_kraken2.sbatch` | default | Kraken2 on `${sample}_microbe_R1.fastq` / `_R2.fastq`; writes minimizer reports |
| 03 | `03_bracken_merge.sbatch` | default | Bracken species estimates, species BIOM, minimizer stats, total-read lookup, Kraken2 QC |
| 05 | `05_downstream.sbatch` -> `05_aggregate_filter.R`, `06_visualize.R` | default | aggregate/filter Kraken2/Bracken species calls, optional decontam, plots, manifest |

---

## Key knobs in `config.sh`

- Kraken2/Bracken paths:
  `KRAKEN2_DB=/oak/stanford/groups/emoding/analysis/brian/tools/microbiome/db/k2_pluspf_20260226`
  and `KRAKEN2_SIF=/oak/stanford/groups/emoding/analysis/brian/tools/microbiome/containers/kraken2_bracken_2.17.1_3.1.sif`.
- Kraken2 defaults: `KRAKEN2_CONFIDENCE=0`,
  `KRAKEN2_MIN_HIT_GROUPS=2`, `REPORT_MINIMIZER_DATA=true`.
- Bracken: `READ_LEN=100`, `BRACKEN_THRESH=10`, `BRACKEN_LEVELS="S"`.
- Abundance/prevalence: `MIN_READS=5`, `MIN_RPM=0.5`,
  `MIN_PREVALENCE_FRAC=0.01`.
- Minimizer filter: `MIN_DISTINCT_MINIMIZERS=5`,
  `MIN_MINIMIZER_COVERAGE=0.01`. Bracken-redistributed species without a raw
  Kraken2 minimizer row still fail open, and the summary reports how often that
  happened.
- Contaminants: `CONTAMINANT_BLOCKLIST=contaminants_blocklist.txt`; optional
  `CONTROLS_LIST` enables `decontam` prevalence filtering.

Step 05 defines `high_confidence` as:

```text
pass_abundance & pass_minimizer & pass_prevalence &
!is_host & !is_blocklisted
```

---

## Outputs (`05_results/`)

- `cross_validated_species.csv` - all Kraken2/Bracken sample-species calls with
  reads, RPM, minimizer metrics, prevalence, blocklist/host flags, and
  `high_confidence`.
- `high_confidence_taxa.csv` - the filtered final calls.
- `cohort_high_confidence_summary.csv` - per-species prevalence and RPM summary.
- `pipeline_summary.txt` - one-page run summary with read funnel, taxon funnel,
  minimizer fail-open counts, and top high-confidence species.
- `QC_read_accounting.csv` - per-sample read accounting plus Kraken2 classified
  and unclassified counts when available.
- `filter_funnel_summary.csv` - counts of sample-taxon calls surviving each gate.
- `rejected_calls.csv.gz` - every call that failed at least one gate, with
  `fail_reason`.
- `run_manifest.txt` - parameter values, tool/DB versions, DB mtimes, and git
  commit state.
- `decontam_contaminants.csv` - only if controls were supplied and `decontam`
  ran.

Figures (`06_figures/`): read-accounting QC, filter funnel, and cohort
high-confidence landscape.

---

## Dependencies

Default steps require the Kraken2/Bracken Apptainer/Singularity container,
the PlusPF Kraken2/Bracken database with `database100mers.kmer_distrib`,
`kraken-biom`, and R 4.4.2 with `readr`, `dplyr`, `tidyr`, `tibble`, `stringr`,
`purrr`, `ggplot2`, `scales`, `phyloseq`, and `biomformat`.

Optional manual step 01 reruns still require STAR, samtools, BBMap tools, bowtie2,
and the host/rRNA/T2T references configured in `config.sh`.
