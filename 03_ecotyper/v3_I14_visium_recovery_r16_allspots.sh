#!/bin/bash
# =============================================================================
# v3_I14_visium_recovery_r16_allspots.sh — EcoTyper Visium recovery of the v3 model on ALL round-16 CytoSPACE spots,
# with the round-16 rule only (author 2026-11-10: "Visium (Fig 5b-e, S5b-d): use v3 rule, re run is possible").
# Same as v3_I4_visium_recovery_r16.sh except for the spot set:
#   v3_I4  : spots of the v1 EcoTyper input (6,091; 81 already lost by the v1 rule) minus round-16 Cycling > 0.5
#   v3_I14 : every spot of the round-16 CytoSPACE run (6,172) minus round-16 Cycling > 0.5
#
# Job 1: v3_I14_visium_fractions_r16_allspots.R (13-type spot fractions + 10x spot folders of the kept spots).
# Job 2: EcoTyper_recovery_visium.R for the 4 tumours, one after the other, with:
#   - Discovery = 2026-10-28_cHL_v2_r16_c0.933_m4_e6_final (the v3 model, as v3_I4)
#   - spot folders 2026-11-10_Visium_Input_r16_allspots/<id>; background "NULL"; 14 threads
#   - spot fractions 2026-11-10_Visium_fractions_r16_allspots/<id>_fractions_r16_allspots.txt
#   - NEW recovery names 2026-11-10_r16_c0.933_m4_e6_visium_r16cyto_allspots_<id> (earlier recoveries are not touched)
#
# Run on Sherlock (from cHL_v2/script/v3):  bash v3_I14_visium_recovery_r16_allspots.sh
# =============================================================================
set -euo pipefail

# Step 1. Settings
MODEL=r16_c0.933_m4_e6
DISCOVERY_NAME=2026-10-28_cHL_v2_${MODEL}_final
ET=${PROJECT_ROOT}/Ecotyper
ECOTYPER_DIR=$ET/ecotyper-master
VISIUM_DIR=$ET/visium/2026-11-10_Visium_Input_r16_allspots         # written by job 1
FRACTIONS_DIR=$ET/visium/2026-11-10_Visium_fractions_r16_allspots  # written by job 1
V2_DIR=${PROJECT_ROOT}/cHL_v2
SCRIPT_DIR=$V2_DIR/script/v3
WORK_DIR=$V2_DIR/output/r16/stepI/visium_recovery_allspots
LOG_DIR=$V2_DIR/output/r16/stepI/logs
TUMOURS="tumor1 tumor2 tumor3 tumor4"
COMMON_SLURM_OPTS="--partition=emoding,normal --ntasks=1 --mail-type=FAIL --mail-user=${MAIL_USER}"
mkdir -p $WORK_DIR $LOG_DIR

# Step 2. Safety checks: the model exists; nothing of this recovery exists yet
if [ ! -d $ECOTYPER_DIR/EcoTyper/$DISCOVERY_NAME ]; then
  echo "discovery $DISCOVERY_NAME is missing" >&2
  exit 1
fi
for tumour in $TUMOURS; do
  if [ -e $ECOTYPER_DIR/2026_H2/2026-11-10_${MODEL}_visium_r16cyto_allspots_${tumour}_Output ]; then
    echo "recovery of $tumour exists; refusing to overwrite" >&2
    exit 1
  fi
done
if [ -e $VISIUM_DIR ] || [ -e $FRACTIONS_DIR ]; then
  echo "an input folder of job 1 exists; refusing to overwrite" >&2
  exit 1
fi

# Step 3. Job 1: spot fractions and spot folders (32 GB: the spot count tables are read whole)
job_fractions=$(sbatch --parsable $COMMON_SLURM_OPTS --cpus-per-task=1 --mem=32G \
  --time=0-01:00:00 --job-name=I14_visium_fractions --output=$LOG_DIR/%x_%j.out <<EOF
#!/bin/bash
module load R/4.2.0
Rscript $SCRIPT_DIR/v3_I14_visium_fractions_r16_allspots.R
EOF
)
echo "I14 fractions: $job_fractions"

# Step 4. Job 2: one YAML per tumour, run one after the other (14 CPUs x 14 GB, 1 day; as v3_I4)
visium_lines=""
for tumour in $TUMOURS; do
  recovery_name=2026-11-10_${MODEL}_visium_r16cyto_allspots_${tumour}
  yaml_file=$WORK_DIR/${recovery_name}.yml
  cat > $yaml_file <<EOF_Y
default :
  Input :
    Discovery dataset name : "${DISCOVERY_NAME}"
    Recovery dataset name : "${recovery_name}"
    Input Visium directory : "${VISIUM_DIR}/${tumour}"
    Recovery cell type fractions : "${FRACTIONS_DIR}/${tumour}_fractions_r16_allspots.txt"
    Background cell type : "NULL"
    CIBERSORTx username : "not_used_user_fractions"
    CIBERSORTx token : "not_used_user_fractions"
  Output :
    Output folder : "${ECOTYPER_DIR}/2026_H2/${recovery_name}_Output"
  Pipeline settings :
    Number of threads : 14
    CIBERSORTx fractions Singularity path : "${ET}/hires.sif"
EOF_Y
  visium_lines="$visium_lines
echo \"=== $tumour \$(date)\"
Rscript EcoTyper_recovery_visium.R -c $yaml_file"
done
job_visium=$(sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$job_fractions --cpus-per-task=14 \
  --mem-per-cpu=14G --time=1-00:00:00 --job-name=I14_visium_recovery --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail
module load system imagemagick/7.0.7-2
module load R/4.2.0
cd $ECOTYPER_DIR
$visium_lines
EOF_J
)
echo "I14 Visium recovery: $job_visium"
