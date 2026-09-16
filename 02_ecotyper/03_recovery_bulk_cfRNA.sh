#!/bin/bash
#SBATCH --job-name=Ecotyper_recovery_cfRNA
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=emoding
# =============================================================================
# 03_recovery_bulk_cfRNA.sh — recover cell states/ecotypes in plasma cfRNA (13 cHL,
# 20 healthy volunteers, 24 NLPHL). Input = cfRNA gene counts normalised as
# TMMwsp -> CPM -> log2 -> RUVg (edgeR::calcNormFactors(method="TMMwsp"), RUVSeq::RUVg),
# genes x samples. No annotation file is needed.
# =============================================================================
module load R/4.2.0
DATASETS="${DATASETS:-/path/to/datasets}"

Rscript EcoTyper_recovery_bulk.R -d cHL_discovery \
  -m ${DATASETS}/cHL_cfRNA_TMMwsp_cpm_log2_ruvg.txt \
  -o cfRNA_recovery
