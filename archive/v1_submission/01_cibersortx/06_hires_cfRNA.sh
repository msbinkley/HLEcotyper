#!/bin/bash
# =============================================================================
# 06_hires_cfRNA.sh
# CIBERSORTx HiRes (S-mode) on plasma cell-free RNA (cfRNA) TPM.
# Mixture = cfRNA TPM (STAR + RSEM), protein-coding genes, columns renormalised to 1e6.
# Outputs feed (i) the cfRNA HRS-fraction comparison and (ii) EcoTyper recovery in cfRNA
# (02_ecotyper/03_recovery_bulk_cfRNA.sh). The same command is used for the LOD95
# spike-in mixtures (04_cfRNA_LOD95/02_cibersortx_hires_spikein.sh).
#
# Usage:  bash 06_hires_cfRNA.sh
# =============================================================================

CIBERSORT_ROOT="${CIBERSORT_ROOT:-/oak/stanford/groups/emoding/analysis/brian/Cibersort}"
INPUT_DIR="${INPUT_DIR:-$CIBERSORT_ROOT/bulk_input}"
HIRES_SIF="${HIRES_SIF:-$CIBERSORT_ROOT/hires.sif}"
PARTITION="${SLURM_PARTITION:-emoding}"
: "${CIBERSORTX_USER:?}"; : "${CIBERSORTX_TOKEN:?}"

mixture="/src/data/cfRNA_coding_symbol_TPM.txt"
sigmatrix="/src/data/cHL_signature_matrix_400_800.txt"
refsample="/src/data/cHL_scRNA_reference_4000.txt"
genelist="/src/data/coding_gene_list.txt"
outdir="$CIBERSORT_ROOT/cHL_cfRNA_400_800/"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_cfRNA
#SBATCH --time=0-02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem-per-cpu=8G
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
