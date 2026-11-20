#!/bin/bash
# =============================================================================
# 04_hires_bulk_tumor.sh
# Deconvolute the 187 FFPE bulk tumor RNA-seq profiles with the cHL signature matrix
# (CIBERSORTx HiRes, S-mode batch correction against the single-cell reference).
#
# This ONE job produces both outputs used downstream:
#   CIBERSORTxGEP_NA_Fractions-Adjusted.txt     -> 16 cell-type fractions per tumor (Fig. 1b)
#   CIBERSORTxHiRes_NA_<celltype>_Window*.txt   -> imputed cell-type-specific expression
#                                                  (input to EcoTyper cell-state discovery)
#
# Inputs (assumed to already exist; see top-level README):
#   mixture    = bulk TPM matrix, genes x samples, protein-coding genes, columns sum to 1e6
#   sigmatrix  = resources/cHL_signature_matrix_400_800.txt
#   refsample  = the single-cell reference used to build the signature matrix
#   genelist   = resources/coding_gene_list.txt  (restricts HiRes imputation to coding genes)
#
# Usage:  bash 04_hires_bulk_tumor.sh
# =============================================================================

CIBERSORT_ROOT="${CIBERSORT_ROOT:-/oak/stanford/groups/emoding/analysis/brian/Cibersort}"
INPUT_DIR="${INPUT_DIR:-$CIBERSORT_ROOT/bulk_input}"
HIRES_SIF="${HIRES_SIF:-$CIBERSORT_ROOT/hires.sif}"
PARTITION="${SLURM_PARTITION:-emoding}"
: "${CIBERSORTX_USER:?}"; : "${CIBERSORTX_TOKEN:?}"

mixture="/src/data/cHL_bulk_TPM.txt"
sigmatrix="/src/data/cHL_signature_matrix_400_800.txt"
refsample="/src/data/cHL_scRNA_reference_4000.txt"
genelist="/src/data/coding_gene_list.txt"
outdir="$CIBERSORT_ROOT/cHL_hires_400_800/"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_hires
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
