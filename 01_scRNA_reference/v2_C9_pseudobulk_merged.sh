#!/bin/bash
# =============================================================================
# v2_C9_pseudobulk_merged.sh
#
# PURPOSE : Build the "Merged" pseudo-bulk set (v2_C9_pseudobulk_merged.R), then deconvolve every mixture with the merged
#           signature (CIBERSORTx Fractions, no batch correction, as v1 / v2_D1 job B).
# RESULTS : output/stepD/pseudobulk_merged/Merged_<type>_pseudo/CIBERSORTx_Results.txt (next to the other five sets).
# RUN     : on Sherlock, bash v2_C9_pseudobulk_merged.sh   (submits two Slurm jobs; the second waits for the first)
# =============================================================================
set -euo pipefail

# Paths on OAK
W=${PROJECT_ROOT}/cHL_v2     # project root
S=$W/script                                              # scripts
D=$W/output/stepD                                        # step-D (CIBERSORTx) outputs and logs
IN=$D/input                                              # mixtures + signature matrices, mounted as /src/data in the container
SIF=${PROJECT_ROOT}/Cibersort/fractions.sif   # CIBERSORTx Fractions container

# Job 1: build the Merged pseudo-bulk mixtures in R (2 CPUs, 64 GB, 4 h)
j1=$(sbatch --parsable \
  --partition=emoding,normal \
  --ntasks=1 \
  --cpus-per-task=2 \
  --mem=64G \
  --time=0-04:00:00 \
  --mail-type=FAIL \
  --mail-user=${MAIL_USER} \
  --job-name=C9pb_merged \
  --output=$D/logs/%x_%j.out \
  --wrap="module load R/4.2.0; Rscript $S/v2_C9_pseudobulk_merged.R")

# Job 2 (starts only if job 1 succeeds; 4 CPUs x 8 GB, 6 h): CIBERSORTx Fractions on every Merged mixture.
# The job script is the text between <<EOF_J and EOF_J. $W, $S, $D, $IN, $SIF are filled in now; the escaped \$ parts
# are evaluated when the job runs. For each mixture file it: hard-links the file into $IN, makes the output folder
# $D/pseudobulk_merged/Merged_<type>_pseudo, and runs CIBERSORTxFractions with signature sigmtx_merged.txt
# (credentials read at run time by v2_D_csx_credentials.sh; csx_mask hides them in the log).
j2=$(sbatch --parsable \
  --dependency=afterok:$j1 \
  --partition=emoding,normal \
  --ntasks=1 \
  --cpus-per-task=4 \
  --mem-per-cpu=8G \
  --time=0-06:00:00 \
  --mail-type=FAIL \
  --mail-user=${MAIL_USER} \
  --job-name=C9csx_merged \
  --output=$D/logs/%x_%j.out <<EOF_J
#!/bin/bash
set -eo pipefail; source $S/v2_D_csx_credentials.sh
for f in $W/output/stepC/pseudobulk/pseudobulk_Merged/pseudobulk_Merged_*_cpm.txt; do
  ln -f \$f $IN/
  n=\$(basename \$f _cpm.txt); n=\${n#pseudobulk_}; out=$D/pseudobulk_merged/\${n}_pseudo; mkdir -p \$out; echo "== \$n"
  singularity exec -c -B $IN:/src/data -B \$out:/src/outdir $SIF /src/CIBERSORTxFractions \\
    --username "\$CSX_USER" --token "\$CSX_TOKEN" --sigmatrix /src/data/sigmtx_merged.txt --mixture /src/data/\$(basename \$f) --verbose TRUE 2>&1 | csx_mask
done
EOF_J
)
echo "pseudo-bulk job $j1 -> CIBERSORTx job $j2"
