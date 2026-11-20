#!/bin/bash
# v3_J1_markers_r16.sh — Slurm job for v3_J1_markers_r16.R (round-16 FindAllMarkers + dot-plot statistics for Fig S1c)
# Run (from OAK cHL_v2/script/v3):  sbatch v3_J1_markers_r16.sh
#SBATCH --job-name=J1_markers_r16
#SBATCH --partition=emoding,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=${PROJECT_ROOT}/cHL_v2/output/r16/stepJ/logs/%x_%j.out
set -euo pipefail
# Step 1. Same R module as the other v2/v3 Sherlock R steps
module load R/4.2.0
# Step 2. Run
Rscript --vanilla ${PROJECT_ROOT}/cHL_v2/script/v3/v3_J1_markers_r16.R
