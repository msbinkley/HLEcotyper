#!/bin/bash
# =============================================================================
# v3_I13b_cibersortx_lod_separate.sh — the two separate S-mode runs of the Fig 6d LOD (author 2026-11-05: "run s mode on
# background group and spike in group seperately"). Mixtures: v3_I13a.
# Recipe as v3_I6b: the ORIGINAL round-16 signature as input + round-16 raw-count reference, S-mode. Each run makes its own
#   adjusted signature from its own samples; no adjusted signature from any other run is used.
# Job script and credential handling copied from v3_I12b unchanged (token read at run time by v2_D_csx_credentials.sh,
#   never printed; output piped through csx_mask; outputs checked and redacted at exit).
# Resources: CIBERSORTx used about one core and at most 22 GB in earlier runs.
#   Background, 94 samples: about 2 to 3 h (146 samples took 2.7 h); limit 12 h.
#   Spiked, 1,504 samples: about 20 h (1,385 samples took 18.9 h); limit 2 days, the most the normal partition allows.
# Run on Sherlock (from cHL_v2/script/v3), after v3_I13a:  bash v3_I13b_cibersortx_lod_separate.sh
# Outputs: OAK cHL_v2/output/r16/stepI/lod_separate_background/ and lod_separate_spiked/CIBERSORTx_Adjusted.txt
# =============================================================================
set -euo pipefail

# Step 1. Folders and files
V2_DIR=${PROJECT_ROOT}/cHL_v2
SCRIPT_DIR=$V2_DIR/script                                   # holds v2_D_csx_credentials.sh
R16_INPUT=$V2_DIR/output/r16/stepD/input                    # original sigmtx_merged.txt and ref_merged.txt of round 16
STEP_DIR=$V2_DIR/output/r16/stepI
INPUT_DIR=$STEP_DIR/lod_separate_input                      # v3_I13a output; bound to /src/data inside the container
LOG_DIR=$STEP_DIR/logs
CIBERSORTX_SIF=${PROJECT_ROOT}/Cibersort/fractions.sif

# Step 2. Round-16 signature and reference, hard-linked next to the mixtures (no copies)
ln -f $R16_INPUT/sigmtx_merged.txt $INPUT_DIR/sigmtx_merged.txt
ln -f $R16_INPUT/ref_merged.txt $INPUT_DIR/ref_merged.txt

# Step 3. One S-mode job per group. Never overwrite an earlier run; the mixture must exist
for group in background spiked; do
  mixture=2026-11-05_lod_$group.txt
  output_dir=$STEP_DIR/lod_separate_$group
  if [ $group = background ]; then
    time_limit=12:00:00
  else
    time_limit=2-00:00:00
  fi
  if [ -e $output_dir ]; then
    echo "$output_dir exists; refusing to overwrite" >&2
    exit 1
  fi
  if [ ! -s $INPUT_DIR/$mixture ]; then
    echo "$INPUT_DIR/$mixture missing; run v3_I13a first" >&2
    exit 1
  fi
  mkdir -p $output_dir $LOG_DIR
  group_job=$(sbatch --parsable --partition=emoding,normal --ntasks=1 --mail-type=FAIL,END --mail-user=${MAIL_USER} \
    --cpus-per-task=6 --mem=32G --time=$time_limit \
    --job-name=I13_lod_$group --output=$LOG_DIR/%x_%j.out <<EOF
#!/bin/bash
set -eo pipefail
source $SCRIPT_DIR/v2_D_csx_credentials.sh
redact () {
  # count output files holding the token (only the count is printed), then replace it ('|| true': no match is fine)
  token_files=\$( (grep -rlIF "\$CSX_TOKEN" $output_dir 2>/dev/null || true) | wc -l)
  echo "output files holding the token: \$token_files"
  (grep -rlIF "\$CSX_TOKEN" $output_dir 2>/dev/null || true) | while read -r token_file; do
    sed -i "s|\$CSX_TOKEN|<TOKEN>|g" "\$token_file"
  done
}
trap redact EXIT
singularity exec -c -B $INPUT_DIR:/src/data -B $output_dir:/src/outdir $CIBERSORTX_SIF /src/CIBERSORTxFractions \\
  --username "\$CSX_USER" --token "\$CSX_TOKEN" --sigmatrix /src/data/sigmtx_merged.txt \\
  --mixture /src/data/$mixture --refsample /src/data/ref_merged.txt --rmbatchSmode TRUE --verbose TRUE 2>&1 | csx_mask
EOF
  )
  echo "I13 S-mode $group: $group_job"
done
