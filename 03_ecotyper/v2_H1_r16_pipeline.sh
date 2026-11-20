#!/bin/bash
# =============================================================================
# v2_H1_r16_pipeline.sh — ROUND 16: reference -> CIBERSORTx -> EcoTyper discovery -> 2-ecotype search with Steidl
# recovery, submitted as ONE chain of Slurm jobs (author 2026-10-26, before going to sleep: "run cibersort and ecotyper
# and identify a two ecotype model that is significant in steidl cohort").
#
# Every step repeats the round-14 recipe; only the input cells (round 16, 110,043 cells) and the output folders change.
# All round-16 files go to output/r16/ and to new EcoTyper names (2026-10-26_cHL_v2_r16_*), so nothing of round 14
# (the current standard) is overwritten.
#
#   J1 reference      v2_C4_reference_pseudobulk.R (merged reference, 4,000 genes, no downsampling, batch filter off,
#                     all sources), round-16 object + labels; SKIP_PSEUDOBULK=1 (the Fig 2c benchmark mixtures are not
#                     needed for EcoTyper / Steidl; they can be made later)                       after the round-16 job
#   J2 signature      CIBERSORTx signature matrix (v1 settings: --G.min 400 --G.max 800, filter TRUE), undotted   after J1
#   J3 bulk fractions CIBERSORTx S-mode on the 187 tumours (refsample = the round-16 reference); pooled to
#                     13 types (v2_D_pool13.R: B_postGC, NK, Cycling dropped) and to 12 types with cDC + pDC = DC
#                     (v2_D17_pool12_dc.R)                                                                   after J2
#   J4 Steidl S-mode  CIBERSORTx S-mode on the Steidl array with the round-16 signature + reference (corrected array,
#                     as v2_E6)                                                                              after J2
#   J5 / J6 discovery EcoTyper discovery from scratch, 13 types (main) and 12 types (DC): 100 NMF restarts, max 6
#                     states, cutoff 0.95, min 6 states, step 7 skipped (as v2_D10 / v2_D18)                 after J3
#   J7 / J8 configs   v2_D11_r14_configs.R on each discovery (cutoff 0.900-1.000 x max states 4-6)     after J5 / J6
#   J9+  search       v2_H2_r16_search_chunk.sh: min states 6 and 5, each pass split into NPARTS parallel chunks;
#                     every 2-ecotype model is recovered in Steidl (raw + corrected)        after J7 / J8 and J4
# The Cox / log-rank tests and the choice of model are done afterwards (v2_H3 on the Mac).
#
# Why min 6 and 5 only: in round 14, 2-ecotype models came from min 6 (main) and min 6 / 5 (DC); min 4 and 3 gave none.
#
# Credentials: read at run time by v2_D_csx_credentials.sh, never printed; CIBERSORTx output piped through csx_mask;
#   EcoTyper discovery YAMLs get the token only inside the job and are redacted on exit (as v2_D10).
#
# How to run (on Sherlock, from cHL_v2/script):  WAIT_JOB=<round-16 job ID> bash v2_H1_r16_pipeline.sh
# =============================================================================
set -euo pipefail

# Step 1. The job to wait for, folders and names
WAIT_JOB=${WAIT_JOB:?set WAIT_JOB (job ID of v2_C22_round16_recluster_and_check.R)}
NPARTS=4                                          # parallel chunks per search pass (round 14: one job per pass, 4.7 h)
SEARCH_MINS="6 5"
V2_DIR=${PROJECT_ROOT}/cHL_v2
SCRIPT_DIR=$V2_DIR/script
R16_DIR=$V2_DIR/output/r16
STEP_C_DIR=$R16_DIR/stepC
STEP_D_DIR=$R16_DIR/stepD
INPUT_DIR=$STEP_D_DIR/input                        # bound to /src/data inside the CIBERSORTx container
SIG_DIR=$STEP_D_DIR/sigmtx_merged
FRACTIONS_DIR=$STEP_D_DIR/fractions_merged
STEIDL_SMODE_DIR=$R16_DIR/stepE/steidl_smode
LOG_DIR=$R16_DIR/logs
ECOTYPER_DIR=${PROJECT_ROOT}/Ecotyper/ecotyper-master
CIBERSORTX_SIF=${PROJECT_ROOT}/Cibersort/fractions.sif
OBJ=cHL_merged_harmony_v2_r16.rds
LABELS=$V2_DIR/output/stepC/inputs/2026-10-26_round16_labels.tsv.gz
BULK_FILE=$V2_DIR/ref/bulk_TPM_v47_coding_187_symbol.tsv
STEIDL_RAW=$V2_DIR/output/stepE/input/2026-10-22_Steidl_microarray_v47_linear.txt
ANNOTATION_FILE=$V2_DIR/output/stepD/input/ecotyper_annotation_187.txt
MAIN_DISCOVERY=2026-10-26_cHL_v2_r16_merged_discovery
DC_DISCOVERY=2026-10-26_cHL_v2_r16_DC_discovery

# Step 2. Safety checks: never write into an existing round-16 folder or discovery; fixed inputs exist
if [ -e $R16_DIR ]; then
  echo "$R16_DIR exists; refusing to mix runs" >&2
  exit 1
fi
for discovery in $MAIN_DISCOVERY $DC_DISCOVERY; do
  if [ -e $ECOTYPER_DIR/EcoTyper/$discovery ]; then
    echo "discovery $discovery exists; refusing to mix runs" >&2
    exit 1
  fi
done
for input_file in $BULK_FILE $STEIDL_RAW $ANNOTATION_FILE $SCRIPT_DIR/v2_H2_r16_search_chunk.sh; do
  if [ ! -s $input_file ]; then
    echo "missing $input_file" >&2
    exit 1
  fi
done
mkdir -p $STEP_C_DIR $INPUT_DIR $SIG_DIR $FRACTIONS_DIR $STEIDL_SMODE_DIR $LOG_DIR \
         $STEP_D_DIR/ecotyper_r16_merged $STEP_D_DIR/ecotyper_r16_DC

# Step 3. Settings shared by the jobs
COMMON_SLURM_OPTS="--partition=emoding,normal --ntasks=1 --mail-type=FAIL --mail-user=${MAIL_USER}"
JOB_HEADER='set -eo pipefail; source '$SCRIPT_DIR'/v2_D_csx_credentials.sh'

# Step 4. J1: reference (v2_C4 through v2_run.sh; 8 CPUs, 128 GB as v2_run.sh)
job_reference=$(sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$WAIT_JOB \
  --cpus-per-task=8 --mem=128G --time=0-12:00:00 --job-name=H1_reference --output=$LOG_DIR/%x_%j.out \
  --export=ALL,OBJ=$OBJ,LABELS=$LABELS,OUT_SUBDIR=output/r16/stepC,SKIP_PSEUDOBULK=1 \
  $SCRIPT_DIR/v2_run.sh v2_C4_reference_pseudobulk.R)
echo "J1 reference: $job_reference"

# Step 5. J2: signature matrix (as v2_D1 job A); the inputs are hard-linked into the container's input folder first
job_signature=$(sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$job_reference \
  --cpus-per-task=10 --mem-per-cpu=20G --time=0-10:00:00 --job-name=H1_signature --output=$LOG_DIR/%x_%j.out <<EOF_A
#!/bin/bash
$JOB_HEADER
ln -f $STEP_C_DIR/reference/cHL_v2_scRNA_reference_4000_merged.txt $INPUT_DIR/ref_merged.txt
ln -f $BULK_FILE $INPUT_DIR/bulk_TPM_v47_coding_187_symbol.tsv
ln -f $STEIDL_RAW $INPUT_DIR/steidl_v47_linear.txt
echo "reference cells: \$(head -1 $INPUT_DIR/ref_merged.txt | awk -F'\t' '{print NF - 1}')"
singularity exec -c -B $INPUT_DIR:/src/data -B $SIG_DIR:/src/outdir $CIBERSORTX_SIF /src/CIBERSORTxFractions \\
  --username "\$CSX_USER" --token "\$CSX_TOKEN" --refsample /src/data/ref_merged.txt \\
  --single_cell TRUE --fraction 0 --filter TRUE --G.min 400 --G.max 800 2>&1 | csx_mask
raw=\$(ls $SIG_DIR/CIBERSORTx_*inferred_phenoclasses*.bm.K999.txt | head -1)
python3 $SCRIPT_DIR/v2_D_undot_sigmatrix.py $INPUT_DIR/ref_merged.txt "\$raw" $SIG_DIR/cHL_v2_sigmtx_400_800_merged.txt
cp $SIG_DIR/cHL_v2_sigmtx_400_800_merged.txt $INPUT_DIR/sigmtx_merged.txt
echo "signature genes: \$(( \$(wc -l < $INPUT_DIR/sigmtx_merged.txt) - 1 ))"
EOF_A
)
echo "J2 signature: $job_signature"

# Step 6. J3: bulk S-mode fractions (as v2_D1 job C; round 14 took 5 h), then the 13-type and 12-type (DC) pooling
job_bulk=$(sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$job_signature \
  --cpus-per-task=8 --mem-per-cpu=8G --time=1-00:00:00 --job-name=H1_bulk --output=$LOG_DIR/%x_%j.out <<EOF_C
#!/bin/bash
$JOB_HEADER
singularity exec -c -B $INPUT_DIR:/src/data -B $FRACTIONS_DIR:/src/outdir $CIBERSORTX_SIF /src/CIBERSORTxFractions \\
  --username "\$CSX_USER" --token "\$CSX_TOKEN" --sigmatrix /src/data/sigmtx_merged.txt \\
  --mixture /src/data/bulk_TPM_v47_coding_187_symbol.tsv --refsample /src/data/ref_merged.txt --rmbatchSmode TRUE --verbose TRUE 2>&1 | csx_mask
module load R/4.2.0
Rscript $SCRIPT_DIR/v2_D_pool13.R $FRACTIONS_DIR/CIBERSORTx_Adjusted.txt $FRACTIONS_DIR/cHL_v2_fractions13_merged.txt
Rscript $SCRIPT_DIR/v2_D17_pool12_dc.R $FRACTIONS_DIR/CIBERSORTx_Adjusted.txt $FRACTIONS_DIR/cHL_v2_fractions12_DC_merged.txt
EOF_C
)
echo "J3 bulk fractions: $job_bulk"

# Step 7. J4: Steidl S-mode (as v2_E6 job 1; round 14 took 3.6 h)
job_steidl=$(sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$job_signature \
  --cpus-per-task=8 --mem-per-cpu=8G --time=1-00:00:00 --job-name=H1_steidl_smode --output=$LOG_DIR/%x_%j.out <<EOF_S
#!/bin/bash
$JOB_HEADER
singularity exec -c -B $INPUT_DIR:/src/data -B $STEIDL_SMODE_DIR:/src/outdir $CIBERSORTX_SIF /src/CIBERSORTxFractions \\
  --username "\$CSX_USER" --token "\$CSX_TOKEN" --sigmatrix /src/data/sigmtx_merged.txt \\
  --mixture /src/data/steidl_v47_linear.txt --refsample /src/data/ref_merged.txt --rmbatchSmode TRUE --verbose TRUE 2>&1 | csx_mask
EOF_S
)
echo "J4 Steidl S-mode: $job_steidl"

# Step 8. J5 / J6: EcoTyper discoveries from scratch (as v2_D10 / v2_D18; round 14 took 3 h each)
#   $1 = discovery name, $2 = fractions file, $3 = run folder (template YAML), $4 = output folder, $5 = job name
#   The template YAML has placeholders for the credentials; the job fills them in and redacts the token on exit.
submit_discovery () {
  local discovery_name=$1
  local fractions_file=$2
  local run_dir=$3
  local output_folder=$4
  local job_name=$5
  local template_yaml=$run_dir/${discovery_name}.template.yaml
  cat > $template_yaml <<EOF_Y
default :
  Input :
    Discovery dataset name : "${discovery_name}"
    Expression matrix : "${BULK_FILE}"
    Cell type fractions : "${fractions_file}"
    Expression type : "RNA-seq"
    Annotation file : "${ANNOTATION_FILE}"
    Annotation file column to scale by : "NULL"
    Annotation file column(s) to plot : [ "EBV", "H1H2", "Age_Bin", "Cohort", "Gender", "Subtype_Code"]
    CIBERSORTx username : "@CSX_USER@"
    CIBERSORTx token : "@CSX_TOKEN@"
  Output :
    Output folder : "${output_folder}"
  Pipeline settings :
    Pipeline steps to skip : [7]
    Number of threads : 24
    Number of NMF restarts : 100
    Maximum number of states per cell type : 6
    Cophenetic coefficient cutoff : 0.95
    CIBERSORTx fractions Singularity path : "${ECOTYPER_DIR}/../fractions.sif"
    CIBERSORTx hires Singularity path : "${ECOTYPER_DIR}/../hires.sif"
    Minimum number of states in ecotypes : 6
EOF_Y
  sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$job_bulk \
    --cpus-per-task=24 --mem-per-cpu=4G --time=2-00:00:00 --job-name=$job_name --output=$LOG_DIR/%x_%j.out <<EOF_J
#!/bin/bash
$JOB_HEADER
if [ ! -s $fractions_file ]; then
  echo "missing $fractions_file" >&2
  exit 1
fi
YAML_FILE=$ECOTYPER_DIR/${discovery_name}.yaml
sed -e "s|@CSX_USER@|\$CSX_USER|" -e "s|@CSX_TOKEN@|\$CSX_TOKEN|" $template_yaml > \$YAML_FILE
redact () {
  grep -rlIF "\$CSX_TOKEN" \$YAML_FILE \\
       $ECOTYPER_DIR/${output_folder} \\
       $ECOTYPER_DIR/EcoTyper/${discovery_name} \\
       $ECOTYPER_DIR/CIBERSORTx/hires/${discovery_name} \\
       $ECOTYPER_DIR/CIBERSORTx/fractions/discovery/${discovery_name} \\
       $ECOTYPER_DIR/datasets/discovery/${discovery_name} 2>/dev/null |
  while read -r token_file; do
    sed -i "s|\$CSX_TOKEN|<TOKEN>|g" "\$token_file"
  done
}
trap redact EXIT
module load R/4.2.0
cd $ECOTYPER_DIR
Rscript EcoTyper_discovery_bulk.R -c \$YAML_FILE 2>&1 | csx_mask
EOF_J
}
job_discovery_main=$(submit_discovery $MAIN_DISCOVERY $FRACTIONS_DIR/cHL_v2_fractions13_merged.txt \
  $STEP_D_DIR/ecotyper_r16_merged 2026_H2/2026-10-26_cHL_v2_r16_merged_r100_95_6_6 H1_discovery_main)
job_discovery_dc=$(submit_discovery $DC_DISCOVERY $FRACTIONS_DIR/cHL_v2_fractions12_DC_merged.txt \
  $STEP_D_DIR/ecotyper_r16_DC 2026_H2/2026-10-26_cHL_v2_r16_DC_r100_95_6_6 H1_discovery_dc)
echo "J5 discovery main: $job_discovery_main; J6 discovery DC: $job_discovery_dc"

# Step 9. J7 / J8: search configurations (v2_D11 with the discovery name and output file as arguments)
submit_configs () {
  local discovery_name=$1
  local run_dir=$2
  local after_job=$3
  local job_name=$4
  sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$after_job \
    --cpus-per-task=1 --mem=8G --time=0-01:00:00 --job-name=$job_name --output=$LOG_DIR/%x_%j.out <<EOF_K
#!/bin/bash
set -eo pipefail
module load R/4.2.0
Rscript $SCRIPT_DIR/v2_D11_r14_configs.R $discovery_name $run_dir/r16_configs.tsv
EOF_K
}
job_configs_main=$(submit_configs $MAIN_DISCOVERY $STEP_D_DIR/ecotyper_r16_merged $job_discovery_main H1_configs_main)
job_configs_dc=$(submit_configs $DC_DISCOVERY $STEP_D_DIR/ecotyper_r16_DC $job_discovery_dc H1_configs_dc)
echo "J7 configs main: $job_configs_main; J8 configs DC: $job_configs_dc"

# Step 10. J9+: search chunks (8 CPUs x 6 GB as v2_D12). They wait for the configurations AND the Steidl S-mode
#   (afterany for S-mode: if it failed, the chunks still run and recover on the raw array only)
for set_name in main dc; do
  if [ "$set_name" = "main" ]; then
    configs_job=$job_configs_main
  else
    configs_job=$job_configs_dc
  fi
  for min_states in $SEARCH_MINS; do
    for part in $(seq 1 $NPARTS); do
      chunk_job=$(sbatch --parsable $COMMON_SLURM_OPTS --dependency=afterok:$configs_job,afterany:$job_steidl \
        --cpus-per-task=8 --mem-per-cpu=6G --time=1-00:00:00 \
        --job-name=H2_${set_name}_min${min_states}_p${part} --output=$LOG_DIR/%x_%j.out \
        --export=ALL,SET=$set_name,MIN=$min_states,PART=$part,NPARTS=$NPARTS \
        $SCRIPT_DIR/v2_H2_r16_search_chunk.sh)
      echo "search $set_name min $min_states part $part: $chunk_job"
    done
  done
done
