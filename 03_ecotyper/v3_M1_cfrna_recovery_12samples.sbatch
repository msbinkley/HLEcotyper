#!/bin/bash
#SBATCH --job-name=M1_cfRNA_recovery_12
#SBATCH --partition=emoding,normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-01:00:00
#SBATCH --output=M1_cfRNA_recovery_12_%j.out
# =============================================================================
# v3_M1_cfrna_recovery_12samples.sbatch — EcoTyper bulk recovery of the final model (r16_c0.933_m4_e6) in plasma
# cfRNA of the 12 patients with cHL used for cell-state recovery (Fig. 6f, left).
#
# Input: linear TPM matrix of these 12 plasma samples (17,865 genes x 12 samples; genes in rows, first column "Gene").
# Command, as for the other bulk recoveries (no -a / -c, default threads):
#   EcoTyper_recovery_bulk.R -d <discovery> -m <matrix> -o <output>
# This run is the cfRNA recovery of the paper. The cfRNA step inside v2_H4_model_chain.sh was an earlier run on an
# earlier version of the matrix and is not used.
# EcoTyper scales every gene across the samples of the matrix, so the result of a sample depends on the other
# samples of the matrix: recover exactly the cohort that is analysed.
#
# How to run (on the cluster, from this folder):  sbatch v3_M1_cfrna_recovery_12samples.sbatch
# Randomness: none.
# =============================================================================
set -euo pipefail

# Step 1. Folders and names
ET=${PROJECT_ROOT}/Ecotyper
ECOTYPER_DIR=$ET/ecotyper-master
DISCOVERY_NAME=2026-10-28_cHL_v2_r16_c0.933_m4_e6_final
CFRNA_MATRIX=$ET/datasets/2026-11-24_cHL_cfRNA_TPM_12samples.txt     # linear TPM, 17,865 genes x 12 samples
OUTPUT_DIR=$ECOTYPER_DIR/2026_H2/2026-11-24_r16_c0.933_m4_e6_cfRNA_12samples
COHORT=$(basename $CFRNA_MATRIX .txt)                                # EcoTyper names the recovery after the matrix file
INTERNAL=$ECOTYPER_DIR/EcoTyper/$DISCOVERY_NAME/Custom

# Step 2. Stop if an output exists already (never overwrite)
if [ -e $OUTPUT_DIR ] || [ -e $INTERNAL/Cell_States/recovery/$COHORT ]; then
  echo "an output of this job exists already; nothing done" >&2
  exit 1
fi

# Step 3. Check the matrix: 12 samples
n_samples=$(head -1 $CFRNA_MATRIX | tr '\t' '\n' | tail -n +2 | wc -l)
echo "matrix: $n_samples samples; $(wc -l < $CFRNA_MATRIX) lines"
if [ "$n_samples" != "12" ]; then
  echo "the matrix does not hold the 12 expected samples" >&2
  exit 1
fi

# Step 4. The recovery. EcoTyper stops with an error when no sample reaches an ecotype; Fig. 6f needs only the
#         cell-state files, so that stop is allowed and handled in step 5.
module load R/4.2.0
cd $ECOTYPER_DIR
echo "=== cfRNA, 12 samples $(date)"
if ! Rscript EcoTyper_recovery_bulk.R -d $DISCOVERY_NAME -m $CFRNA_MATRIX -o $OUTPUT_DIR; then
  echo "EcoTyper_recovery_bulk.R stopped before the end"
fi

# Step 5. If EcoTyper did not copy its results (stop at the ecotype step), copy the cell-state files of every cell
#         type at its chosen number of states, exactly as EcoTyper_recovery_bulk.R does (as v2_H4_bulk_finish.sh).
#         In the run of the paper EcoTyper finished normally and this step did nothing.
FINAL_OUTPUT=$OUTPUT_DIR/$COHORT
if [ -s $FINAL_OUTPUT/Ecotypes/ecotype_assignment.txt ]; then
  echo "complete: $FINAL_OUTPUT"
else
  if [ ! -d $INTERNAL/Cell_States/recovery/$COHORT ]; then
    echo "no cell-state recovery was written; the run failed" >&2
    exit 1
  fi
  mkdir -p $FINAL_OUTPUT/Ecotypes
  tail -n +2 $INTERNAL/Analysis/rank_selection/rank_data.txt | while read -r cell_type n_states rest; do
    cell_type=${cell_type//\"/}                                      # rank_data.txt may quote the cell type
    state_dir=$INTERNAL/Cell_States/recovery/$COHORT/$cell_type/$n_states
    mkdir -p $FINAL_OUTPUT/$cell_type
    for state_file in state_abundances.txt state_assignment.txt state_assignment_heatmap.pdf heatmap_data.txt heatmap_top_ann.txt; do
      if [ -e $state_dir/$state_file ]; then
        cp -n $state_dir/$state_file $FINAL_OUTPUT/$cell_type/
      else
        echo "missing $state_dir/$state_file"
      fi
    done
  done
  for ecotype_file in initial_ecotype_assignment.txt heatmap_all_samples.pdf; do
    if [ -e $INTERNAL/Ecotypes/recovery/$COHORT/$ecotype_file ]; then
      cp -n $INTERNAL/Ecotypes/recovery/$COHORT/$ecotype_file $FINAL_OUTPUT/Ecotypes/
    fi
  done
  echo "EcoTyper assigned no sample of $COHORT to an ecotype ($(date +%F)); cell-state files copied by v3_M1" \
    > $FINAL_OUTPUT/Ecotypes/NO_SAMPLES_ASSIGNED.txt
  echo "cell states copied: $(ls -d $FINAL_OUTPUT/*/ | wc -l) folders in $FINAL_OUTPUT"
fi
echo "=== done $(date)"
