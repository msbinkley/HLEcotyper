#!/bin/bash
#SBATCH --job-name=Ecotyper_recovery_EPIC
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=emoding
# =============================================================================
# 04_recovery_bulk_EPICseq.sh — recover cell states/ecotypes from plasma EPIC-seq
# inferred expression (Alizadeh lab, 2.6-Mb TSS panel, 1,676 genes; public data).
# Input = EPIC-seq inferred expression on a LINEAR scale, genes x samples.
# =============================================================================
module load R/4.2.0
DATASETS="${DATASETS:-/path/to/datasets}"

Rscript EcoTyper_recovery_bulk.R -d cHL_discovery \
  -m ${DATASETS}/cHL_epic_linear_mixture.txt \
  -o EPIC_recovery
