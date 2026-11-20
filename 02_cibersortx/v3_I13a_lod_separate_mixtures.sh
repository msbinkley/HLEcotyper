#!/bin/bash
# v3_I13a_lod_separate_mixtures.sh — Slurm job for v3_I13a_lod_separate_mixtures.py (about 3.5 million noise draws, 28 million mixed values)
# Run (from OAK cHL_v2/script/v3):  sbatch v3_I13a_lod_separate_mixtures.sh
#SBATCH --job-name=I13a_lod_mixtures
#SBATCH --partition=emoding,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=${PROJECT_ROOT}/cHL_v2/output/r16/stepI/logs/%x_%j.out
set -euo pipefail
python3 ${PROJECT_ROOT}/cHL_v2/script/v3/v3_I13a_lod_separate_mixtures.py
