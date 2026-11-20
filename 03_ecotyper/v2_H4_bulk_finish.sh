#!/bin/bash
# =============================================================================
# v2_H4_bulk_finish.sh — finish the cfRNA and EPIC-seq recoveries of one hit model when its B job
# (v2_H4_model_chain.sh, "H4_<model>_bulk") stopped early.
#
# Why: EcoTyper_recovery_bulk.R first recovers the cell states of every cell type, then assigns samples to ecotypes, and
# only then copies the results to the output folder. When no sample of a cohort reaches an ecotype it stops with
# "No samples were assigned to ecotypes!" (or, for EPIC-seq, an empty-heatmap error) and copies nothing; because the B job
# runs with set -e, any recovery after it never starts. On 2026-10-28 this happened for:
#   r16_c0.991_m4_e6    cfRNA ok; EPIC-seq stopped at ecotype assignment
#   r17_c0.979_m5_e6    cfRNA stopped at ecotype assignment; EPIC-seq never ran
#   r17DC_c0.936_m6_e6  cfRNA stopped at ecotype assignment; EPIC-seq never ran
# Fig 6f reads only the per-cell-type state recoveries (heatmap_data.txt, heatmap_top_ann.txt), which exist.
#
# What it does, for cfRNA and then EPIC-seq:
#   - output complete (Ecotypes/ecotype_assignment.txt present)        -> nothing
#   - EcoTyper's own cell-state recovery for the cohort missing        -> run EcoTyper_recovery_bulk.R (a stop at the
#                                                                        ecotype step is allowed and reported)
#   - output still not complete                                         -> copy the cell-state files exactly as
#     EcoTyper_recovery_bulk.R would (same files, same rank key), copy initial_ecotype_assignment.txt and
#     heatmap_all_samples.pdf to Ecotypes/, and write Ecotypes/NO_SAMPLES_ASSIGNED.txt
# Steidl is not touched (it finished for every model).
#
# How to run (on Sherlock, from cHL_v2/script):
#   bash v2_H4_bulk_finish.sh <model> [<job ID to wait for (afterany)>]
# =============================================================================
set -euo pipefail

# Step 1. Model and the job to wait for (a still-pending B job of the same model)
MODEL=${1:?model name}
WAIT_FOR_JOB=${2:-}

# Step 2. Folders and cohort matrices (same as v2_H4_model_chain.sh)
ET=${PROJECT_ROOT}/Ecotyper
ECOTYPER_DIR=$ET/ecotyper-master
V2_DIR=${PROJECT_ROOT}/cHL_v2
LOG_DIR=$V2_DIR/output/hits/logs
DISCOVERY_NAME=2026-10-28_cHL_v2_${MODEL}_final
RECOVERY_BASE=$ECOTYPER_DIR/2026_H2/2026-10-28_${MODEL}
CFRNA_MATRIX=$ET/datasets/2025-08-01_cHL_cfRNA_TPM.txt
EPIC_MATRIX=$ET/datasets/2026-03-10_cHL_epic_linear_mixture.txt
COMMON_SLURM_OPTS="--partition=emoding,normal --ntasks=1 --mail-type=FAIL --mail-user=${MAIL_USER}"

# Step 3. The model must exist
if [ ! -s $ECOTYPER_DIR/EcoTyper/$DISCOVERY_NAME/Custom/Analysis/rank_selection/rank_data.txt ]; then
  echo "model $DISCOVERY_NAME not found" >&2
  exit 1
fi

# Step 4. Wait for the model's B job if it is still pending or running (afterany: run whether it passes or fails)
dependency=""
if [ -n "$WAIT_FOR_JOB" ]; then
  dependency="--dependency=afterany:$WAIT_FOR_JOB"
fi

# Step 5. The job (2 CPUs x 8 GB, 2 h: each recovery took about 2 min with 8 CPUs)
job_finish=$(sbatch --parsable $COMMON_SLURM_OPTS $dependency --cpus-per-task=2 --mem-per-cpu=8G --time=0-02:00:00 \
  --job-name=H4_${MODEL}_bulkfinish --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
set -uo pipefail
module load R/4.2.0
cd $ECOTYPER_DIR
internal=EcoTyper/$DISCOVERY_NAME/Custom
rank_key=\$internal/Analysis/rank_selection/rank_data.txt

for pair in cfRNA:$CFRNA_MATRIX EPIC:$EPIC_MATRIX; do
  cohort_label=\${pair%%:*}
  matrix_file=\${pair#*:}
  cohort=\$(basename \$matrix_file .txt)                     # EcoTyper names the recovery after the matrix file
  final_output=${RECOVERY_BASE}_\$cohort_label/\$cohort
  echo "=== \$cohort_label (\$cohort) \$(date)"

  # Step 5a. Complete already: nothing to do
  if [ -s \$final_output/Ecotypes/ecotype_assignment.txt ]; then
    echo "complete: \$final_output"
    continue
  fi

  # Step 5b. EcoTyper's cell-state recovery missing (the cohort never ran): run it; a stop at the ecotype step is allowed
  if [ ! -d \$internal/Cell_States/recovery/\$cohort ]; then
    echo "running EcoTyper_recovery_bulk.R"
    if ! Rscript EcoTyper_recovery_bulk.R -d $DISCOVERY_NAME -m \$matrix_file -o ${RECOVERY_BASE}_\$cohort_label; then
      echo "EcoTyper_recovery_bulk.R stopped (expected when no sample reaches an ecotype)"
    fi
    if [ -s \$final_output/Ecotypes/ecotype_assignment.txt ]; then
      echo "complete: \$final_output"
      continue
    fi
  fi

  # Step 5c. Copy the cell-state recovery of every cell type at its chosen rank, as EcoTyper_recovery_bulk.R does
  mkdir -p \$final_output/Ecotypes
  tail -n +2 \$rank_key | while read -r cell_type n_states rest; do
    cell_type=\${cell_type//\"/}                               # rank_data.txt may quote the cell type
    state_dir=\$internal/Cell_States/recovery/\$cohort/\$cell_type/\$n_states
    mkdir -p \$final_output/\$cell_type
    for state_file in state_abundances.txt state_assignment.txt state_assignment_heatmap.pdf heatmap_data.txt heatmap_top_ann.txt; do
      if [ -e \$state_dir/\$state_file ]; then
        cp -f \$state_dir/\$state_file \$final_output/\$cell_type/
      else
        echo "missing \$state_dir/\$state_file"
      fi
    done
  done

  # Step 5d. The ecotype step's partial output, and a note that no sample was assigned to an ecotype
  for ecotype_file in initial_ecotype_assignment.txt heatmap_all_samples.pdf; do
    if [ -e \$internal/Ecotypes/recovery/\$cohort/\$ecotype_file ]; then
      cp -f \$internal/Ecotypes/recovery/\$cohort/\$ecotype_file \$final_output/Ecotypes/
    fi
  done
  echo "EcoTyper assigned no sample of \$cohort to an ecotype with model $MODEL (\$(date +%F)); cell-state files copied by v2_H4_bulk_finish.sh" \
    > \$final_output/Ecotypes/NO_SAMPLES_ASSIGNED.txt
  echo "cell states copied: \$(ls -d \$final_output/*/ | wc -l) folders in \$final_output"
done
EOF_J
)
echo "bulk finish for $MODEL: $job_finish"
