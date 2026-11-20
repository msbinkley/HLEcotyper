#!/bin/bash
# =============================================================================
# v3_I2_cytospace_r16.sh — CytoSPACE on the 4 Visium tumours with the ROUND-16 scRNA reference (author 2026-10-30:
# "do both"; "use all 110,043 r16 cells").
#
# Same command as v1 (Cytospace/2025-07-12_cHL_cytospace.sh): CytoSPACE v1.1.0, default settings (seed 1, Seurat
# label-transfer cell-type fractions, mean 5 cells per spot, scRNA transcripts downsampled), solver lap_CSPR, the same
# v1 Visium inputs (ST_data_<id>.csv, Coordinates_<id>.csv). Only the scRNA inputs change (v3_I1).
# Memory: v1 used 8 CPUs x 8 GB for 44,769 cells; round 16 has 110,043 cells, so 8 CPUs and 200 GB.
#
# Runs (Sherlock, from cHL_v2/script/v3):  bash v3_I2_cytospace_r16.sh [<job id of v3_I1>]
# Outputs: OAK Cytospace/output_r16/<id>_results/ (v1 outputs in Cytospace/output/ are not touched)
# =============================================================================
set -euo pipefail

# Step 1. Inputs and outputs
CYTOSPACE_DIR=${PROJECT_ROOT}/Cytospace
INPUT_DIR=${PROJECT_ROOT}/cHL_v2/output/r16/stepI/cytospace_input
SCRNA_DATA=$INPUT_DIR/2026-10-30_r16_scRNA_data.txt
SCRNA_LABELS=$INPUT_DIR/2026-10-30_r16_cell_type_labels.txt
OUTPUT_ROOT=$CYTOSPACE_DIR/output_r16
LOG_DIR=${PROJECT_ROOT}/cHL_v2/output/r16/stepI/logs
TUMOURS="tumor1 tumor2 tumor3 tumor4"
WAIT_JOB=${1:-}
mkdir -p $OUTPUT_ROOT $LOG_DIR

# Step 2. Wait for the input job when one is given
dependency=""
if [ -n "$WAIT_JOB" ]; then
  dependency="--dependency=afterok:$WAIT_JOB"
fi

# Step 3. One job per tumour (stop if CytoSPACE results exist: never overwrite)
#   (2026-10-30: the first submission failed at "conda activate" before CytoSPACE started and left empty <id>_results/
#    folders, so the check is on the results subfolder)
for tumour in $TUMOURS; do
  outdir=$OUTPUT_ROOT/${tumour}_results
  if [ -e $outdir/cytospace_results ]; then
    echo "$outdir exists; refusing to overwrite" >&2
    exit 1
  fi
  sbatch --parsable $dependency --job-name=I2_cyto_$tumour --partition=emoding,normal --ntasks=1 --cpus-per-task=8 \
    --mem=200G --time=1-00:00:00 --mail-type=FAIL --mail-user=${MAIL_USER} --output=$LOG_DIR/%x_%j.out <<EOF
#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate $CYTOSPACE_DIR/.conda_envs/cytospace_v1.1.0     # full path (the env is not in the default env list)
mkdir -p $outdir
cd $outdir
cytospace \\
   --scRNA-path $SCRNA_DATA \\
   --cell-type-path $SCRNA_LABELS \\
   --st-path $CYTOSPACE_DIR/cHL_Visium_Input/ST_data_${tumour}.csv \\
   --coordinates-path $CYTOSPACE_DIR/cHL_Visium_Input/Coordinates_${tumour}.csv \\
   --solver-method lap_CSPR
EOF
done
