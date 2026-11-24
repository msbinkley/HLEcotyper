#!/bin/bash
# v3_I10a_mean_sorted_hrs16.sh — Slurm job for v3_I10a_mean_sorted_hrs16.py (reads four GTFs; ~10 min, as v3_I9a)
# Run (from OAK cHL_v2/script/v3):  sbatch v3_I10a_mean_sorted_hrs16.sh
#SBATCH --job-name=I10a_mean_hrs16
#SBATCH --partition=emoding,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/cHL_v2/output/r16/stepI/logs/%x_%j.out
set -euo pipefail
python3 ${PROJECT_ROOT}/cHL_v2/script/v3/v3_I10a_mean_sorted_hrs16.py
