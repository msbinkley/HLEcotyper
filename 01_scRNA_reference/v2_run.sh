#!/bin/bash
#SBATCH --job-name=cHL_v2
#SBATCH --partition=emoding
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
# =============================================================================
# v2_run.sh — runs one step script of the v2 re-analysis on Sherlock (Slurm settings above: 8 CPUs, 128 GB, 24 h).
#
# How to run (from W/script, W = OAK cHL_v2):
#   sbatch --output="../logs/$(date +%F)_<step>_%j.out" v2_run.sh <script>.R
# =============================================================================

# Step 1. Load R 4.2.0 and go to the script folder
module load R/4.2.0
cd ${PROJECT_ROOT}/cHL_v2/script

# Step 2. Run the R script given as the first argument
Rscript "$1"
