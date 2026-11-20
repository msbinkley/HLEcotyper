#!/bin/bash
# =============================================================================
# v2_G6_recovery_exvivo.sh — EcoTyper scRNA recovery of the FINAL round-14 model in the ex vivo Flex samples (Fig 8e);
# author 2026-10-25. Same settings as v1 (2026_H2/2026-10-12_ecotyper_recovery_exvivo.sh):
#   EcoTyper_recovery_scRNA.R -d <discovery> -m <count> -a <annotation> -o <output> -c Sample -t 8 -z TRUE
# One job: inputs (v2_G6_exvivo_recovery_inputs.R), then the recovery.
#
# How to run (on Sherlock):  bash v2_G6_recovery_exvivo.sh
# =============================================================================
set -euo pipefail

# Step 1. Folders and names
ET=${PROJECT_ROOT}/Ecotyper
ECOTYPER_DIR=$ET/ecotyper-master
DISCOVERY_NAME=2026-10-24_cHL_v2_r14_final
V2_DIR=${PROJECT_ROOT}/cHL_v2
COUNT_FILE=$ET/datasets/2026-10-25_cHL_r14_exvivo_count.txt
ANNOT_FILE=$ET/datasets/2026-10-25_cHL_r14_exvivo_annotation.txt
OUTPUT_DIR=$ECOTYPER_DIR/2026_H2/2026-10-25_cHL_r14_exvivo_recovery

# Step 2. Never write into an existing output
if [ -e $OUTPUT_DIR ]; then
  echo "$OUTPUT_DIR exists; refusing to overwrite" >&2
  exit 1
fi

# Step 3. Submit (8 CPUs x 16 GB, 2 days)
sbatch --parsable --partition=emoding,normal --ntasks=1 \
  --cpus-per-task=8 --mem-per-cpu=16G --time=2-00:00:00 --mail-type=END,FAIL --mail-user=${MAIL_USER} \
  --job-name=G6recov_exvivo --output=$V2_DIR/output/stepE/logs/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail
module load R/4.2.0
Rscript $V2_DIR/script/v2_G6_exvivo_recovery_inputs.R
cd $ECOTYPER_DIR
Rscript EcoTyper_recovery_scRNA.R -d $DISCOVERY_NAME -m $COUNT_FILE -a $ANNOT_FILE -o $OUTPUT_DIR -c Sample -t 8 -z TRUE
EOF_J
