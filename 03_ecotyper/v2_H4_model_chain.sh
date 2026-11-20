#!/bin/bash
# =============================================================================
# v2_H4_model_chain.sh — build ONE hit model as a permanent EcoTyper model and run every recovery its figure set needs
# (author 2026-10-28: "for the 8 hits, create a folder for each hits under [figures_v2], and generate entire fig 2 to
# fig 8 (except for fig7) and supp figures for each hit ... always use corrected steidl"; Fig 3e with -s 2500).
#
# The same steps that made the round-14 final model and its recoveries, with the hit's settings:
#   A  final model      as v2_D15: a dedicated copy of the round's discovery, steps 4-8 with the hit's cutoff / max / min
#                       states; checks: chosen ranks = the search vector, Ecotypes files identical to the search run
#   B  bulk recoveries  (after A) Steidl on the round's CIBERSORTx S-mode (corrected) array, -a / -c as the search chunk;
#                       cfRNA and EPIC-seq as v2_G2
#   C  Visium           (after A) as v2_G3, the 4 tumours with the round-14 spot fractions (same pools); DC model:
#                       DC = DC_Conventional + DC_Plasmacytoid, written in the job
#   D  ex vivo          (after A) as v2_G6, the round-14 ex vivo input (same cells and labels); DC model: DC_Conventional -> DC
#   E  merged scRNA     (after A, and after the round's input job when given) as v2_G4 but with -s 2500 (author
#                       2026-10-28) and 4 cell types at a time (-t 4, so the memory fits; the full run of round 14 ran out
#                       of memory with 24)
# Round 14 (model r14_c0.9465_m4_e6) runs only E, on its existing final model (2026-10-24_cHL_v2_r14_final) and inputs.
#
# Names (all new; nothing existing is overwritten):
#   discovery  EcoTyper/2026-10-28_cHL_v2_<model>_final
#   outputs    2026_H2/2026-10-28_cHL_v2_<model>_final/ (model), 2026_H2/2026-10-28_<model>_<recovery>/ (recoveries)
#
# How to run (on Sherlock, from cHL_v2/script):
#   bash v2_H4_model_chain.sh <model> [<job ID of the round's scRNA input job>]
# Models: r14_c0.9465_m4_e6, r16_c0.9275_m4_e6, r16_c0.933_m4_e6, r16_c0.991_m4_e6, r17_c0.979_m5_e6, r17_c0.994_m4_e6,
#         r17_c0.9675_m6_e6, r17_c0.971_m5_e5, r17DC_c0.936_m6_e6
# =============================================================================
set -euo pipefail

# Step 1. The model's settings (from the search tables figures_v2/EcoTyper_search/*_search_ALL.tsv: tag, vec)
MODEL=${1:?model name}
SCRNA_INPUT_JOB=${2:-}
case $MODEL in
  r14_c0.9465_m4_e6)  ROUND=r14; SET=main; CUTOFF=0.9465; MAX_STATES=4; MIN_STATES=6; RANKS=3-3-2-3-3-4-2-2-3-3-4-3 ;;
  r16_c0.9275_m4_e6)  ROUND=r16; SET=main; CUTOFF=0.9275; MAX_STATES=4; MIN_STATES=6; RANKS=4-3-4-4-4-2-3-2-3-3-2-2-3 ;;
  r16_c0.933_m4_e6)   ROUND=r16; SET=main; CUTOFF=0.933;  MAX_STATES=4; MIN_STATES=6; RANKS=4-3-3-4-4-2-3-2-3-3-2-2-3 ;;
  r16_c0.991_m4_e6)   ROUND=r16; SET=main; CUTOFF=0.991;  MAX_STATES=4; MIN_STATES=6; RANKS=3-4-2-3-3-2-3-2-4-3-2-2-3 ;;
  r17_c0.979_m5_e6)   ROUND=r17; SET=main; CUTOFF=0.979;  MAX_STATES=5; MIN_STATES=6; RANKS=4-4-2-3-2-3-4-4-4-4-4-3 ;;
  r17_c0.994_m4_e6)   ROUND=r17; SET=main; CUTOFF=0.994;  MAX_STATES=4; MIN_STATES=6; RANKS=3-3-2-3-3-3-2-3-2-3-3-3 ;;
  r17_c0.9675_m6_e6)  ROUND=r17; SET=main; CUTOFF=0.9675; MAX_STATES=6; MIN_STATES=6; RANKS=5-5-2-4-5-4-4-4-4-4-5-4 ;;
  r17_c0.971_m5_e5)   ROUND=r17; SET=main; CUTOFF=0.971;  MAX_STATES=5; MIN_STATES=5; RANKS=3-5-2-4-2-4-4-4-4-4-5-3 ;;
  r17DC_c0.936_m6_e6) ROUND=r17; SET=dc;   CUTOFF=0.936;  MAX_STATES=6; MIN_STATES=6; RANKS=5-3-2-5-5-5-6-4-4-4-4-4 ;;
  *) echo "unknown model $MODEL" >&2; exit 1 ;;
esac
TAG=c${CUTOFF}_m${MAX_STATES}_e${MIN_STATES}

# Step 2. Folders and input files
ET=${PROJECT_ROOT}/Ecotyper
ECOTYPER_DIR=$ET/ecotyper-master
V2_DIR=${PROJECT_ROOT}/cHL_v2
EXPRESSION_FILE=$V2_DIR/ref/bulk_TPM_v47_coding_187_symbol.tsv
ANNOTATION_FILE=$V2_DIR/output/stepD/input/ecotyper_annotation_187.txt
STEIDL_CLINICAL=$V2_DIR/output/stepE/input/Steidl_clinical.txt
CFRNA_MATRIX=$ET/datasets/2025-08-01_cHL_cfRNA_TPM.txt                      # as v2_G2
EPIC_MATRIX=$ET/datasets/2026-03-10_cHL_epic_linear_mixture.txt             # as v2_G2
V1_VISIUM_DIR=$ET/visium/2025-07-21_Visium_Input_Cyto                       # as v2_G3
VISIUM_FRACTIONS_DIR=$ET/visium/2026-10-24_Visium_fractions_r14             # round-14 pools (B_postGC, NK = NK + ILC)
TUMOURS="tumor1 tumor2 tumor3 tumor4"
EXVIVO_COUNT=$ET/datasets/2026-10-25_cHL_r14_exvivo_count.txt               # as v2_G6
EXVIVO_ANNOT=$ET/datasets/2026-10-25_cHL_r14_exvivo_annotation.txt
LOG_DIR=$V2_DIR/output/hits/logs
WORK_DIR=$V2_DIR/output/hits/$MODEL                                          # YAMLs and the DC inputs of this model
mkdir -p $LOG_DIR $WORK_DIR
COMMON_SLURM_OPTS="--partition=emoding,normal --ntasks=1 --mail-type=FAIL --mail-user=${MAIL_USER}"

if [ "$ROUND" = "r14" ]; then
  DISCOVERY_NAME=2026-10-24_cHL_v2_r14_final                                 # the existing final model (v2_D15)
  SCRNA_COUNT=$ET/datasets/2026-10-25_cHL_r14_scRNA_merged_count.txt         # v2_G4 inputs
  SCRNA_ANNOT=$ET/datasets/2026-10-25_cHL_r14_scRNA_merged_annotation.txt
else
  DISCOVERY_NAME=2026-10-28_cHL_v2_${MODEL}_final
  ROUND_DIR=$V2_DIR/output/$ROUND
  if [ "$SET" = "main" ]; then
    SOURCE_DISCOVERY=2026-10-26_cHL_v2_${ROUND}_merged_discovery
    FRACTIONS_FILE=$ROUND_DIR/stepD/fractions_merged/cHL_v2_fractions13_merged.txt
    SEARCH_RUN=$ECOTYPER_DIR/2026_H2/2026-10-26_cHL_v2_${ROUND}_merged_search/min${MIN_STATES}/$TAG
  else
    SOURCE_DISCOVERY=2026-10-26_cHL_v2_${ROUND}_DC_discovery
    FRACTIONS_FILE=$ROUND_DIR/stepD/fractions_merged/cHL_v2_fractions12_DC_merged.txt
    SEARCH_RUN=$ECOTYPER_DIR/2026_H2/2026-10-26_cHL_v2_${ROUND}_DC_search/min${MIN_STATES}/$TAG
  fi
  STEIDL_CORRECTED=$ROUND_DIR/stepE/steidl_smode/CIBERSORTx_Mixtures_Adjusted.txt
  SCRNA_COUNT=$ET/datasets/2026-10-28_cHL_${ROUND}_scRNA_merged_count.txt     # v2_H4_round_extras.sh J1
  if [ "$SET" = "dc" ]; then
    SCRNA_ANNOT=$ET/datasets/2026-10-28_cHL_${ROUND}_scRNA_merged_annotation_DC.txt
  else
    SCRNA_ANNOT=$ET/datasets/2026-10-28_cHL_${ROUND}_scRNA_merged_annotation.txt
  fi
fi
MODEL_OUTPUT=2026_H2/2026-10-28_cHL_v2_${MODEL}_final                          # relative to the EcoTyper folder
RECOVERY_BASE=$ECOTYPER_DIR/2026_H2/2026-10-28_${MODEL}

# Step 3. Safety checks: nothing of this model exists yet
for existing in $ECOTYPER_DIR/$MODEL_OUTPUT ${RECOVERY_BASE}_steidl_corrected ${RECOVERY_BASE}_cfRNA ${RECOVERY_BASE}_EPIC \
                ${RECOVERY_BASE}_exvivo ${RECOVERY_BASE}_scRNA_s2500; do
  if [ "$ROUND" != "r14" ] || [ "$existing" = "${RECOVERY_BASE}_scRNA_s2500" ]; then
    if [ -e $existing ]; then
      echo "$existing exists; refusing to overwrite" >&2
      exit 1
    fi
  fi
done

# Step 4. A: the final model (not for round 14, which has one)
dependency=""
if [ "$ROUND" != "r14" ]; then
  if [ -e $ECOTYPER_DIR/EcoTyper/$DISCOVERY_NAME ] || [ ! -s $SEARCH_RUN/Ecotypes/ecotypes.txt ]; then
    echo "discovery $DISCOVERY_NAME exists, or the search run $SEARCH_RUN is missing" >&2
    exit 1
  fi
  YAML_FILE=$WORK_DIR/${DISCOVERY_NAME}.yaml
  cat > $YAML_FILE <<EOF_Y
default :
  Input :
    Discovery dataset name : "${DISCOVERY_NAME}"
    Expression matrix : "${EXPRESSION_FILE}"
    Cell type fractions : "${FRACTIONS_FILE}"
    Expression type : "RNA-seq"
    Annotation file : "${ANNOTATION_FILE}"
    Annotation file column to scale by : "NULL"
    Annotation file column(s) to plot : [ "EBV", "H1H2", "Age_Bin", "Cohort", "Gender", "Subtype_Code"]
    CIBERSORTx username : "not_used_steps_1_2_skipped"
    CIBERSORTx token : "not_used_steps_1_2_skipped"
  Output :
    Output folder : "${MODEL_OUTPUT}"
  Pipeline settings :
    Pipeline steps to skip : [1,2,3,7]
    Number of threads : 8
    Number of NMF restarts : 100
    Maximum number of states per cell type : ${MAX_STATES}
    Cophenetic coefficient cutoff : ${CUTOFF}
    CIBERSORTx fractions Singularity path : "${ECOTYPER_DIR}/../fractions.sif"
    CIBERSORTx hires Singularity path : "${ECOTYPER_DIR}/../hires.sif"
    Minimum number of states in ecotypes : ${MIN_STATES}
EOF_Y
  job_model=$(sbatch --parsable $COMMON_SLURM_OPTS --cpus-per-task=8 --mem-per-cpu=6G --time=0-06:00:00 \
    --job-name=H4_${MODEL}_model --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail
module load R/4.2.0
cd $ECOTYPER_DIR
rsync -a --exclude=restarts EcoTyper/$SOURCE_DISCOVERY/ EcoTyper/$DISCOVERY_NAME/
cp -al CIBERSORTx/hires/$SOURCE_DISCOVERY CIBERSORTx/hires/$DISCOVERY_NAME
cp -al CIBERSORTx/fractions/discovery/$SOURCE_DISCOVERY CIBERSORTx/fractions/discovery/$DISCOVERY_NAME
echo "copied discovery to $DISCOVERY_NAME"
Rscript EcoTyper_discovery_bulk.R -c $YAML_FILE
got=\$(tail -n +2 EcoTyper/$DISCOVERY_NAME/Custom/Analysis/rank_selection/rank_data.txt | cut -f2 | paste -sd- -)
if [ "\$got" != "$RANKS" ]; then
  echo "RANK MISMATCH: got \$got" >&2
  exit 1
fi
echo "ranks ok"
for f in ecotypes.txt ecotype_assignment.txt ecotype_abundance.txt; do
  if cmp -s $MODEL_OUTPUT/Ecotypes/\$f $SEARCH_RUN/Ecotypes/\$f; then
    echo "\$f identical to the search run"
  else
    echo "\$f DIFFERS from the search run" >&2
    exit 1
  fi
done
EOF_J
)
  echo "A model: $job_model"
  dependency="--dependency=afterok:$job_model"

  # Step 5. B: Steidl (corrected array), cfRNA, EPIC-seq (8 CPUs x 8 GB, 12 h; as v2_G2 and the search chunk)
  job_bulk=$(sbatch --parsable $COMMON_SLURM_OPTS $dependency --cpus-per-task=8 --mem-per-cpu=8G --time=0-12:00:00 \
    --job-name=H4_${MODEL}_bulk --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail
module load R/4.2.0
cd $ECOTYPER_DIR
echo "=== Steidl corrected \$(date)"
mkdir -p ${RECOVERY_BASE}_steidl_corrected
Rscript EcoTyper_recovery_bulk.R -d $DISCOVERY_NAME -m $STEIDL_CORRECTED -a $STEIDL_CLINICAL -c stage,Age45 -t 8 \
  -o ${RECOVERY_BASE}_steidl_corrected
echo "=== cfRNA \$(date)"
Rscript EcoTyper_recovery_bulk.R -d $DISCOVERY_NAME -m $CFRNA_MATRIX -o ${RECOVERY_BASE}_cfRNA
echo "=== EPIC-seq \$(date)"
Rscript EcoTyper_recovery_bulk.R -d $DISCOVERY_NAME -m $EPIC_MATRIX -o ${RECOVERY_BASE}_EPIC
EOF_J
)
  echo "B bulk recoveries: $job_bulk"

  # Step 6. C: Visium, the 4 tumours one after the other (14 CPUs x 14 GB, 1 day; as v2_G3)
  #   DC model: spot fractions with DC = DC_Conventional + DC_Plasmacytoid (other columns unchanged)
  if [ "$SET" = "dc" ]; then
    fractions_dir=$WORK_DIR/visium_fractions_DC
    mkdir -p $fractions_dir
    for tumour in $TUMOURS; do
      awk -F'\t' -v OFS='\t' '
        NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "DC_Conventional") cdc = i; if ($i == "DC_Plasmacytoid") pdc = i } }
        { line = ""
          for (i = 1; i <= NF; i++) {
            if (i == pdc) continue
            value = $i
            if (i == cdc) value = (NR == 1) ? "DC" : sprintf("%.15g", $cdc + $pdc)    # full precision (awk default is 6 digits)
            line = (line == "") ? value : line OFS value
          }
          print line }' $VISIUM_FRACTIONS_DIR/${tumour}_fractions_r14.txt > $fractions_dir/${tumour}_fractions_DC.txt
    done
    fractions_suffix=_fractions_DC.txt
  else
    fractions_dir=$VISIUM_FRACTIONS_DIR
    fractions_suffix=_fractions_r14.txt
  fi
  visium_lines=""
  for tumour in $TUMOURS; do
    recovery_name=2026-10-28_${MODEL}_visium_${tumour}
    yaml_file=$WORK_DIR/${recovery_name}.yml
    cat > $yaml_file <<EOF_Y
default :
  Input :
    Discovery dataset name : "${DISCOVERY_NAME}"
    Recovery dataset name : "${recovery_name}"
    Input Visium directory : "${V1_VISIUM_DIR}/${tumour}"
    Recovery cell type fractions : "${fractions_dir}/${tumour}${fractions_suffix}"
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
  job_visium=$(sbatch --parsable $COMMON_SLURM_OPTS $dependency --cpus-per-task=14 --mem-per-cpu=14G --time=1-00:00:00 \
    --job-name=H4_${MODEL}_visium --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail
module load system imagemagick/7.0.7-2
module load R/4.2.0
cd $ECOTYPER_DIR
$visium_lines
EOF_J
)
  echo "C Visium: $job_visium"

  # Step 7. D: ex vivo scRNA recovery (8 CPUs x 16 GB, 2 days; as v2_G6). DC model: DC_Conventional -> DC
  if [ "$SET" = "dc" ]; then
    exvivo_annotation=$WORK_DIR/exvivo_annotation_DC.txt
    awk -F'\t' -v OFS='\t' 'NR > 1 && $2 == "DC_Conventional" { $2 = "DC" } { print }' $EXVIVO_ANNOT > $exvivo_annotation
  else
    exvivo_annotation=$EXVIVO_ANNOT
  fi
  job_exvivo=$(sbatch --parsable $COMMON_SLURM_OPTS $dependency --cpus-per-task=8 --mem-per-cpu=16G --time=2-00:00:00 \
    --job-name=H4_${MODEL}_exvivo --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail
module load R/4.2.0
cd $ECOTYPER_DIR
Rscript EcoTyper_recovery_scRNA.R -d $DISCOVERY_NAME -m $EXVIVO_COUNT -a $exvivo_annotation -o ${RECOVERY_BASE}_exvivo -c Sample -t 8 -z TRUE
EOF_J
)
  echo "D ex vivo: $job_exvivo"
fi

# Step 8. E: merged scRNA recovery with -s 2500 (8 CPUs x 24 GB, 2 days; 4 cell types at a time)
if [ -n "$SCRNA_INPUT_JOB" ]; then
  if [ -n "$dependency" ]; then
    dependency="$dependency:$SCRNA_INPUT_JOB"
  else
    dependency="--dependency=afterok:$SCRNA_INPUT_JOB"
  fi
fi
job_scrna=$(sbatch --parsable $COMMON_SLURM_OPTS $dependency --cpus-per-task=8 --mem-per-cpu=24G --time=2-00:00:00 \
  --job-name=H4_${MODEL}_scRNA --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail
module load R/4.2.0
cd $ECOTYPER_DIR
Rscript EcoTyper_recovery_scRNA.R -d $DISCOVERY_NAME -m $SCRNA_COUNT -a $SCRNA_ANNOT -o ${RECOVERY_BASE}_scRNA_s2500 \
  -c Dataset -t 4 -z TRUE -s 2500
EOF_J
)
echo "E scRNA (-s 2500): $job_scrna"
