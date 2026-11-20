#!/bin/bash
# =============================================================================
# v2_H2_r16_search_chunk.sh — one chunk of the ROUND-16 EcoTyper search (author 2026-10-26: "run cibersort and
# ecotyper and identify a two ecotype model that is significant in steidl cohort"). Submitted by
# v2_H1_r16_pipeline.sh as a Slurm job; everything is decided at run time, because the configurations only exist
# after the discovery has finished.
#
# Same search as round 14 (v2_D11 configurations + v2_D12 runs + v2_D19 Steidl recovery), with two differences:
#   1. The configurations of one min-states pass are split over NPARTS jobs (configuration number modulo NPARTS), each
#      in its OWN copy of the discovery (<discovery>_min<M>_p<part>), so they run in parallel. (EcoTyper keeps one
#      internal model per discovery name, so runs inside one copy must be sequential.)
#   2. Right after a configuration gives exactly 2 ecotypes, that model is recovered in the Steidl cohort on BOTH
#      arrays (raw = uncorrected linear array, as Fig 4f; corrected = CIBERSORTx S-mode with the round-16 signature).
#      The internal model of the copy is that model at that moment, so no rerun of the discovery is needed
#      (v2_D19 reran it because its recoveries came after the whole search).
#
# Settings per configuration (as v2_D12): steps 1, 2, 3, 7 skipped (no CIBERSORTx call, no token), 100 NMF restarts
#   (reused from the discovery), the configuration's representative max states and cophenetic cutoff, MIN states.
#
# Environment (set by v2_H1): SET (main | dc), MIN (min states per ecotype), PART (1..NPARTS), NPARTS
# Outputs (EcoTyper folder):
#   2026_H2/2026-10-26_cHL_v2_r16_<set>_search/min<M>/<tag>/              (search run of each configuration)
#   2026_H2/2026-10-26_cHL_v2_r16_<set>_steidl/<tag>/steidl_{raw,corrected}/   (2-ecotype models only)
#   chunk log: output/r16/logs/H2_<set>_min<M>_p<part>_<job>.out (one line per configuration: ranks / ecotypes / n_eco)
# =============================================================================
set -uo pipefail     # no -e: one failing configuration must not stop the others
# (every Rscript inside the loop reads stdin from /dev/null, so it cannot swallow lines of the configuration file)

# Step 1. Which chunk
SET=${SET:?}
MIN=${MIN:?}
PART=${PART:?}
NPARTS=${NPARTS:?}

# Step 2. Folders and names
V2_DIR=${PROJECT_ROOT}/cHL_v2
R16_DIR=$V2_DIR/output/r16
STEP_D_DIR=$R16_DIR/stepD
ECOTYPER_DIR=${PROJECT_ROOT}/Ecotyper/ecotyper-master
EXPRESSION_FILE=$V2_DIR/ref/bulk_TPM_v47_coding_187_symbol.tsv
ANNOTATION_FILE=$V2_DIR/output/stepD/input/ecotyper_annotation_187.txt     # unchanged clinical annotation (read only)
STEIDL_RAW=$V2_DIR/output/stepE/input/2026-10-22_Steidl_microarray_v47_linear.txt
STEIDL_CORRECTED=$R16_DIR/stepE/steidl_smode/CIBERSORTx_Mixtures_Adjusted.txt
STEIDL_CLINICAL=$V2_DIR/output/stepE/input/Steidl_clinical.txt
if [ "$SET" = "main" ]; then
  SOURCE_DISCOVERY=2026-10-26_cHL_v2_r16_merged_discovery
  FRACTIONS_FILE=$STEP_D_DIR/fractions_merged/cHL_v2_fractions13_merged.txt
  RUN_DIR=$STEP_D_DIR/ecotyper_r16_merged
  SEARCH_BASE=2026_H2/2026-10-26_cHL_v2_r16_merged_search/min${MIN}      # relative to the EcoTyper folder
  STEIDL_BASE=2026_H2/2026-10-26_cHL_v2_r16_merged_steidl
elif [ "$SET" = "dc" ]; then
  SOURCE_DISCOVERY=2026-10-26_cHL_v2_r16_DC_discovery
  FRACTIONS_FILE=$STEP_D_DIR/fractions_merged/cHL_v2_fractions12_DC_merged.txt
  RUN_DIR=$STEP_D_DIR/ecotyper_r16_DC
  SEARCH_BASE=2026_H2/2026-10-26_cHL_v2_r16_DC_search/min${MIN}
  STEIDL_BASE=2026_H2/2026-10-26_cHL_v2_r16_DC_steidl
else
  echo "SET must be main or dc" >&2
  exit 1
fi
DISCOVERY_NAME=${SOURCE_DISCOVERY}_min${MIN}_p${PART}                  # this chunk's own copy
RANK_SELECTION_DIR=$ECOTYPER_DIR/EcoTyper/$DISCOVERY_NAME/Custom/Analysis/rank_selection
CONFIG_FILE=$RUN_DIR/r16_configs.tsv
YAML_DIR=$RUN_DIR/search_yaml/min${MIN}
LOG_DIR=$R16_DIR/logs/search_${SET}
mkdir -p $YAML_DIR $LOG_DIR $ECOTYPER_DIR/$SEARCH_BASE $ECOTYPER_DIR/$STEIDL_BASE

# Step 3. Safety checks: configurations exist; Steidl arrays exist (the corrected one only if the S-mode job worked)
if [ ! -s $CONFIG_FILE ]; then
  echo "missing $CONFIG_FILE" >&2
  exit 1
fi
HAVE_CORRECTED=1
if [ ! -s $STEIDL_CORRECTED ]; then
  echo "no corrected Steidl array ($STEIDL_CORRECTED): raw recovery only"
  HAVE_CORRECTED=0
fi

# Step 4. This chunk's copy of the discovery (as v2_D12: without the per-restart NMF folders; the read-only
#         CIBERSORTx fraction and HiRes folders hard-linked under the new name)
module load R/4.2.0
cd $ECOTYPER_DIR
if [ ! -d EcoTyper/$DISCOVERY_NAME ]; then
  rsync -a --exclude=restarts EcoTyper/$SOURCE_DISCOVERY/ EcoTyper/$DISCOVERY_NAME/
  cp -al CIBERSORTx/hires/$SOURCE_DISCOVERY CIBERSORTx/hires/$DISCOVERY_NAME
  cp -al CIBERSORTx/fractions/discovery/$SOURCE_DISCOVERY CIBERSORTx/fractions/discovery/$DISCOVERY_NAME
  echo "copied discovery to $DISCOVERY_NAME"
fi

# Step 5. The configurations of this chunk: number i (1, 2, ...) belongs to part ((i - 1) mod NPARTS) + 1
#   columns of r16_configs.tsv: vec, n_settings, maxes, rep_max, rep_cut, total_states (header line skipped)
config_number=0
while IFS=$'\t' read -r vec n_settings maxes rep_max rep_cut total_states; do
  if [ "$vec" = "vec" ]; then
    continue
  fi
  config_number=$((config_number + 1))
  if [ $(( (config_number - 1) % NPARTS + 1 )) -ne $PART ]; then
    continue
  fi
  tag=c${rep_cut}_m${rep_max}_e${MIN}
  yaml_file=$YAML_DIR/$tag.yaml

  # 5a. YAML (as v2_D12)
  cat > $yaml_file <<EOF_Y
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
    Output folder : "${SEARCH_BASE}/${tag}"
  Pipeline settings :
    Pipeline steps to skip : [1,2,3,7]
    Number of threads : 8
    Number of NMF restarts : 100
    Maximum number of states per cell type : ${rep_max}
    Cophenetic coefficient cutoff : ${rep_cut}
    CIBERSORTx fractions Singularity path : "${ECOTYPER_DIR}/../fractions.sif"
    CIBERSORTx hires Singularity path : "${ECOTYPER_DIR}/../hires.sif"
    Minimum number of states in ecotypes : ${MIN}
EOF_Y

  # 5b. Discovery steps 4-8 -> check the chosen ranks -> count the ecotypes
  echo "=== $tag (config $config_number, expected $vec) $(date)"
  if ! Rscript EcoTyper_discovery_bulk.R -c $yaml_file < /dev/null > $LOG_DIR/min${MIN}_${tag}_discovery.log 2>&1; then
    echo "DISCOVERY FAILED"
    continue
  fi
  got=$(tail -n +2 $RANK_SELECTION_DIR/rank_data.txt | cut -f2 | paste -sd- -)
  if [ "$got" = "$vec" ]; then
    echo "ranks ok"
  else
    echo "RANK MISMATCH: got $got"
  fi
  ecotypes_file=$SEARCH_BASE/$tag/Ecotypes/ecotypes.txt
  if [ ! -f $ecotypes_file ]; then
    echo "no ecotypes"
    continue
  fi
  n_eco=$(tail -n +2 $ecotypes_file | cut -f5 | sort -u | wc -l)
  echo "ecotypes ok; n_eco $n_eco"

  # 5c. Exactly 2 ecotypes: recover this model (the copy's internal model right now) in Steidl, raw and corrected
  #     (EcoTyper creates only the last folder level of its output, so the model folder is made first)
  if [ "$n_eco" -eq 2 ]; then
    mkdir -p $STEIDL_BASE/$tag
    Rscript EcoTyper_recovery_bulk.R -d $DISCOVERY_NAME -m $STEIDL_RAW -a $STEIDL_CLINICAL -c stage,Age45 -t 8 \
      -o $STEIDL_BASE/$tag/steidl_raw < /dev/null > $LOG_DIR/min${MIN}_${tag}_raw.log 2>&1 && echo "steidl raw ok" || echo "STEIDL RAW FAILED"
    if [ $HAVE_CORRECTED -eq 1 ]; then
      Rscript EcoTyper_recovery_bulk.R -d $DISCOVERY_NAME -m $STEIDL_CORRECTED -a $STEIDL_CLINICAL -c stage,Age45 -t 8 \
        -o $STEIDL_BASE/$tag/steidl_corrected < /dev/null > $LOG_DIR/min${MIN}_${tag}_corrected.log 2>&1 && echo "steidl corrected ok" || echo "STEIDL CORRECTED FAILED"
    fi
  fi
done < $CONFIG_FILE
echo "chunk done: $SET min $MIN part $PART of $NPARTS $(date)"
