#!/bin/bash
# =============================================================================
# v3_I6b_cibersortx_cfrna_validation_r16.sh — CIBERSORTx fractions of the CardioGuardian cfRNA WITH the validation
# cohort (NovaSeq B35, CC69-CC97) and the ROUND-16 signature (author 2026-11-02: "rerun deconvolution with the
# validation cohort, and compare cHL vs NLPHL vs other cancers vs healthy").
#
# Mixture: v3_I6a output cfRNA_TPM_141_harmonized_r16.txt (112 earlier + 29 B35 samples, gene names matched to the
#   round-16 signature). All 141 in ONE run, because S-mode batch correction depends on which mixtures are run
#   together (as the 2026-09-25 v1-signature B35 run, Cibersort/2026-09-25_cHL_cfRNA_B35).
# Deconvolution: the same settings as v3_I3 (Fig 6b r16) and the round-16 tumour fractions (v2_H1 J3):
#   CIBERSORTxFractions, round-16 signature (sigmtx_merged.txt), round-16 merged reference (raw counts), S-mode.
#   v3_I3's 112-sample run took 2 h 15 min with 1 CPU busy; 24 h is the limit asked for here.
# Credentials: read at run time by v2_D_csx_credentials.sh, never printed; output piped through csx_mask; after the
#   run every output file is checked for the token (only the count is printed) and redacted.
#
# Run on Sherlock (from cHL_v2/script/v3), after v3_I6a:  bash v3_I6b_cibersortx_cfrna_validation_r16.sh
# Output: OAK cHL_v2/output/r16/stepI/cfrna_validation_fractions/CIBERSORTx_Adjusted.txt
# =============================================================================
set -euo pipefail

# Step 1. Folders and files
V2_DIR=${PROJECT_ROOT}/cHL_v2
SCRIPT_DIR=$V2_DIR/script                                   # holds v2_D_csx_credentials.sh
R16_INPUT=$V2_DIR/output/r16/stepD/input                    # sigmtx_merged.txt, ref_merged.txt of round 16
STEP_DIR=$V2_DIR/output/r16/stepI
INPUT_DIR=$STEP_DIR/cfrna_validation_input                  # v3_I6a output; bound to /src/data inside the container
OUTPUT_DIR=$STEP_DIR/cfrna_validation_fractions
LOG_DIR=$STEP_DIR/logs
CIBERSORTX_SIF=${PROJECT_ROOT}/Cibersort/fractions.sif
MIXTURE=cfRNA_TPM_141_harmonized_r16.txt

# Step 2. Never overwrite an earlier run; the mixture must exist
if [ -e $OUTPUT_DIR ]; then
  echo "$OUTPUT_DIR exists; refusing to overwrite" >&2
  exit 1
fi
if [ ! -s $INPUT_DIR/$MIXTURE ]; then
  echo "$INPUT_DIR/$MIXTURE missing; run v3_I6a first" >&2
  exit 1
fi
mkdir -p $OUTPUT_DIR $LOG_DIR

# Step 3. Round-16 signature and reference, hard-linked next to the mixture (no copies)
ln -f $R16_INPUT/sigmtx_merged.txt $INPUT_DIR/sigmtx_merged.txt
ln -f $R16_INPUT/ref_merged.txt $INPUT_DIR/ref_merged.txt

# Step 4. One job (8 CPUs x 16 GB, as v3_I3; emoding or normal partition)
job_id=$(sbatch --parsable --partition=emoding,normal --ntasks=1 --mail-type=FAIL,END --mail-user=${MAIL_USER} \
  --cpus-per-task=8 --mem-per-cpu=16G --time=24:00:00 \
  --job-name=I6_cfrna_validation --output=$LOG_DIR/%x_%j.out <<EOF
#!/bin/bash
set -eo pipefail
source $SCRIPT_DIR/v2_D_csx_credentials.sh
redact () {
  # count output files holding the token (only the count is printed), then replace it ('|| true': no match is fine)
  token_files=\$( (grep -rlIF "\$CSX_TOKEN" $OUTPUT_DIR 2>/dev/null || true) | wc -l)
  echo "output files holding the token: \$token_files"
  (grep -rlIF "\$CSX_TOKEN" $OUTPUT_DIR 2>/dev/null || true) | while read -r token_file; do
    sed -i "s|\$CSX_TOKEN|<TOKEN>|g" "\$token_file"
  done
}
trap redact EXIT
singularity exec -c -B $INPUT_DIR:/src/data -B $OUTPUT_DIR:/src/outdir $CIBERSORTX_SIF /src/CIBERSORTxFractions \\
  --username "\$CSX_USER" --token "\$CSX_TOKEN" --sigmatrix /src/data/sigmtx_merged.txt \\
  --mixture /src/data/$MIXTURE --refsample /src/data/ref_merged.txt --rmbatchSmode TRUE --verbose TRUE 2>&1 | csx_mask
EOF
)
echo "I6 cfRNA validation fractions: $job_id"
