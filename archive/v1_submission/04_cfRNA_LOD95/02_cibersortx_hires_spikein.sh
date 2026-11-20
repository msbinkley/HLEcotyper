#!/bin/bash
# =============================================================================
# 02_cibersortx_hires_spikein.sh — deconvolute the simulated cfRNA mixtures.
# Identical CIBERSORTx HiRes S-mode call to 01_cibersortx/06_hires_cfRNA.sh; the only
# output used is CIBERSORTxGEP_NA_Fractions-Adjusted.txt (HRS column) -> 03_fit_LOD95.R.
# Usage:  bash 02_cibersortx_hires_spikein.sh
# =============================================================================
CIBERSORT_ROOT="${CIBERSORT_ROOT:-/oak/stanford/groups/emoding/analysis/brian/Cibersort}"
INPUT_DIR="${INPUT_DIR:-$CIBERSORT_ROOT/bulk_input}"
HIRES_SIF="${HIRES_SIF:-$CIBERSORT_ROOT/hires.sif}"
PARTITION="${SLURM_PARTITION:-emoding}"
: "${CIBERSORTX_USER:?}"; : "${CIBERSORTX_TOKEN:?}"

mixture="/src/data/LOD95_HRS_background_spikein_TPM.txt"     # from 01_make_spikein_mixtures.R
sigmatrix="/src/data/cHL_signature_matrix_400_800.txt"
refsample="/src/data/cHL_scRNA_reference_4000.txt"
genelist="/src/data/coding_gene_list.txt"
outdir="$CIBERSORT_ROOT/cfRNA_LOD95_spikein/"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_LOD95
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=14G
#SBATCH --partition=${PARTITION}

mkdir -p ${outdir}
singularity exec -c -B ${INPUT_DIR}:/src/data -B ${outdir}:/src/outdir \\
  ${HIRES_SIF} /src/CIBERSORTxHiRes \\
  --username ${CIBERSORTX_USER} --token ${CIBERSORTX_TOKEN} \\
  --sigmatrix ${sigmatrix} \\
  --mixture ${mixture} \\
  --refsample ${refsample} \\
  --subsetgenes ${genelist} \\
  --rmbatchSmode TRUE
EOF
