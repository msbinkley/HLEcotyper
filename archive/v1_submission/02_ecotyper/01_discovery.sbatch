#!/bin/bash
#SBATCH --job-name=Ecotyper_discovery
#SBATCH --time=02-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=14G
#SBATCH --partition=emoding
# =============================================================================
# 01_discovery.sbatch — EcoTyper cell-state + ecotype DISCOVERY on bulk cHL RNA-seq.
# Run from inside a clone of EcoTyper (https://github.com/digitalcytometry/ecotyper):
#   cd /path/to/ecotyper-master
#   cp /path/to/HLEcotyper/02_ecotyper/discovery_config.yaml .
#   sbatch /path/to/HLEcotyper/02_ecotyper/01_discovery.sbatch
# Results: <Output folder>/<cell type>/{gene_info,state_assignment,state_abundances}.txt,
#          <Output folder>/Ecotypes/{ecotypes,ecotype_assignment,ecotype_abundance}.txt
# The recovery-ready model is written to EcoTyper/<Discovery dataset name>/ inside the clone.
# =============================================================================
module load R/4.2.0
Rscript EcoTyper_discovery_bulk.R -c discovery_config.yaml
