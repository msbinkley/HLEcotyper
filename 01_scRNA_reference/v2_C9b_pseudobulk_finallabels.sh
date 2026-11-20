#!/bin/bash
# =============================================================================
# v2_C9b_pseudobulk_finallabels.sh
#
# PURPOSE : Rebuild the 10X and Merged pseudo-bulk sets with final labels (v2_C9b_pseudobulk_finallabels.R), then run
#           CIBERSORTx Fractions (no batch correction) with all three signatures (merged, lab_only, merged_noS1) for 10X
#           and with the merged signature for Merged. Old mixtures / results are moved (not deleted) to
#           superseded_B_cell_label/ first.
# RUN     : on Sherlock, bash v2_C9b_pseudobulk_finallabels.sh   (submits two Slurm jobs; the second waits for the first)
# =============================================================================
set -euo pipefail

# Paths on OAK
W=${PROJECT_ROOT}/cHL_v2     # project root
S=$W/script                                              # scripts
D=$W/output/stepD                                        # step-D (CIBERSORTx) outputs and logs
IN=$D/input                                              # mixtures + signature matrices, mounted as /src/data in the container
PB=$W/output/stepC/pseudobulk                            # pseudo-bulk mixture sets
SIF=${PROJECT_ROOT}/Cibersort/fractions.sif   # CIBERSORTx Fractions container

# Step 1. Move the old 10X / Merged mixtures and results (built with the "B_cell" label) into superseded_B_cell_label/
mkdir -p $PB/superseded_B_cell_label $IN/superseded_B_cell_label
# 1a. the two mixture-set folders
for s in 10X Merged; do
  [ -d $PB/pseudobulk_$s ] && mv $PB/pseudobulk_$s $PB/superseded_B_cell_label/
done
# 1b. the hard-linked copies of their mixture files in the CIBERSORTx input folder (no error if there are none)
mv $IN/pseudobulk_10X_*_cpm.txt $IN/pseudobulk_Merged_*_cpm.txt $IN/superseded_B_cell_label/ 2>/dev/null || true
# 1c. the CIBERSORTx result folders of each signature
for a in merged lab_only merged_noS1; do
  mkdir -p $D/pseudobulk_$a/superseded_B_cell_label
  for x in $D/pseudobulk_$a/10X_*_pseudo $D/pseudobulk_$a/Merged_*_pseudo; do
    [ -d "$x" ] && mv "$x" $D/pseudobulk_$a/superseded_B_cell_label/
  done
done

# Step 2. Job 1: rebuild the 10X and Merged mixtures in R (2 CPUs, 64 GB, 4 h)
j1=$(sbatch --parsable \
  --partition=emoding,normal \
  --ntasks=1 \
  --cpus-per-task=2 \
  --mem=64G \
  --time=0-04:00:00 \
  --mail-type=FAIL \
  --mail-user=${MAIL_USER} \
  --job-name=C9b_pb \
  --output=$D/logs/%x_%j.out \
  --wrap="module load R/4.2.0; Rscript $S/v2_C9b_pseudobulk_finallabels.R")

# Step 3. Job 2 (starts only if job 1 succeeds; 4 CPUs x 8 GB, 8 h): CIBERSORTx Fractions.
# The job script is the text between <<EOF_J and EOF_J. $W, $S, $D, $IN, $PB, $SIF are filled in now; the escaped \$
# parts are evaluated when the job runs. run <mixture file> <signature arm> hard-links the mixture into $IN, makes
# $D/pseudobulk_<arm>/<set>_<type>_pseudo and runs CIBERSORTxFractions with sigmtx_<arm>.txt (credentials read at run
# time by v2_D_csx_credentials.sh; csx_mask hides them in the log). 10X mixtures: all three arms; Merged: merged only.
j2=$(sbatch --parsable \
  --dependency=afterok:$j1 \
  --partition=emoding,normal \
  --ntasks=1 \
  --cpus-per-task=4 \
  --mem-per-cpu=8G \
  --time=0-08:00:00 \
  --mail-type=FAIL \
  --mail-user=${MAIL_USER} \
  --job-name=C9b_csx \
  --output=$D/logs/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail; source $S/v2_D_csx_credentials.sh
run () { f=\$1; a=\$2; ln -f \$f $IN/; n=\$(basename \$f _cpm.txt); n=\${n#pseudobulk_}; out=$D/pseudobulk_\$a/\${n}_pseudo; mkdir -p \$out; echo "== \$a \$n"
  singularity exec -c -B $IN:/src/data -B \$out:/src/outdir $SIF /src/CIBERSORTxFractions \\
    --username "\$CSX_USER" --token "\$CSX_TOKEN" --sigmatrix /src/data/sigmtx_\$a.txt --mixture /src/data/\$(basename \$f) --verbose TRUE 2>&1 | csx_mask; }
for f in $PB/pseudobulk_10X/pseudobulk_10X_*_cpm.txt; do for a in merged lab_only merged_noS1; do run \$f \$a; done; done
for f in $PB/pseudobulk_Merged/pseudobulk_Merged_*_cpm.txt; do run \$f merged; done
EOF_J
)
echo "pseudo-bulk job $j1 -> CIBERSORTx job $j2"
