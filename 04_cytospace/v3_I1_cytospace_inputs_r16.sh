#!/bin/bash
#SBATCH --job-name=I1_cyto_inputs
#SBATCH --partition=emoding,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=0-06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${MAIL_USER}
# =============================================================================
# v3_I1_cytospace_inputs_r16.sh — runs v3_I1_cytospace_inputs_r16.R (R 4.2.0, as v2_run.sh).
#   sbatch --output=../../output/r16/stepI/logs/%x_%j.out v3_I1_cytospace_inputs_r16.sh   (from cHL_v2/script/v3)
# =============================================================================
module load R/4.2.0
cd ${PROJECT_ROOT}/cHL_v2/script/v3
Rscript v3_I1_cytospace_inputs_r16.R
