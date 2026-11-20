#!/bin/bash
#SBATCH --job-name=Ecotyper_recovery_microarray
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=emoding
# =============================================================================
# 02_recovery_bulk_microarray.sh — recover cHL cell states/ecotypes in the independent
# 130-case RNA-microarray validation cohort (Steidl et al., NEJM 2010; GEO GSE17920).
# Run from the EcoTyper clone that holds the trained model (EcoTyper/cHL_discovery/).
#   -d  Discovery dataset name (from discovery_config.yaml)
#   -m  expression matrix, genes x samples (log2 microarray intensities)
#   -a  sample annotation (ID column + columns given to -c, plotted as heat-map tracks)
#   -o  output folder
# For cross-platform recovery the expression matrix was first batch-corrected with
# CIBERSORTx B-mode against the cHL signature matrix (use CIBERSORTxGEP_NA_Mixture-Adjusted.txt
# from a Fractions run with --rmbatchBmode TRUE as -m).
# =============================================================================
module load R/4.2.0
DATASETS="${DATASETS:-/path/to/datasets}"

Rscript EcoTyper_recovery_bulk.R -d cHL_discovery \
  -m ${DATASETS}/nejm_2010_microarray.txt \
  -a ${DATASETS}/Steidl_annotation.txt \
  -c Age45,Age60 \
  -o Steidl_recovery
