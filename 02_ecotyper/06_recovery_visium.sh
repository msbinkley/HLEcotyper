#!/bin/bash
#SBATCH --job-name=Ecotyper_recovery_visium
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=14G
#SBATCH --partition=emoding
# =============================================================================
# 06_recovery_visium.sh — recover cell states/ecotypes on the four Visium tumors.
# One YAML per tumor (template: visium_recovery_config_template.yml). Two variants were run,
# differing only in the per-spot cell-type fractions supplied:
#   *_Ciber.yml : fractions from CIBERSORTx HiRes on spots (01_cibersortx/05_hires_visium.sh),
#                 renormalised to sum to 1 per spot
#   *_Cyto.yml  : fractions from CytoSPACE (03_cytospace; fractional_abundances_by_spot.csv),
#                 renormalised to sum to 1 per spot
# "Input Visium directory" = Space Ranger output folder for that tumor.
# =============================================================================
module load imagemagick/7.0.7-2
module load R/4.2.0

for t in 18081 42988 42989 47166; do
  Rscript EcoTyper_recovery_visium.R -c visium_recovery_${t}_Ciber.yml
  Rscript EcoTyper_recovery_visium.R -c visium_recovery_${t}_Cyto.yml
done
