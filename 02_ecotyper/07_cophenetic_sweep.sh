#!/bin/bash
#SBATCH --job-name=Ecotyper_sweep
#SBATCH --time=01-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=14G
#SBATCH --partition=emoding
# =============================================================================
# 07_cophenetic_sweep.sh — sensitivity of state/ecotype discovery to the cophenetic
# coefficient cutoff. Re-runs EcoTyper steps 4-8 (rank selection onward) on the existing
# HiRes output for cutoffs 0.90 ... 1.00 and recovers each model in the microarray cohort,
# so the number of states/ecotypes and their prognostic association can be compared.
# The published model uses cutoff 0.90 (discovery_config.yaml).
# Run from the EcoTyper clone after 01_discovery.sbatch has completed once.
# =============================================================================
module load R/4.2.0
DATASETS="${DATASETS:-/path/to/datasets}"

for cutoff in $(seq 0.90 0.0025 1.00); do
  CONFIG="discovery_sweep_${cutoff}.yaml"
  sed -e "s/Cophenetic coefficient cutoff : .*/Cophenetic coefficient cutoff : ${cutoff}/" \
      -e "s/Output folder : .*/Output folder : \"cHL_output_${cutoff}_6_6\"/" \
      -e "s/Pipeline steps to skip : .*/Pipeline steps to skip : [1,2,3,7]/" \
      discovery_config.yaml > "$CONFIG"

  Rscript EcoTyper_discovery_bulk.R -c "$CONFIG"
  Rscript EcoTyper_recovery_bulk.R -d cHL_discovery \
     -m ${DATASETS}/nejm_2010_microarray.txt \
     -a ${DATASETS}/Steidl_annotation.txt \
     -c Age45,Age60 \
     -o "Steidl_recovery_${cutoff}"
  rm -f "$CONFIG"
done
