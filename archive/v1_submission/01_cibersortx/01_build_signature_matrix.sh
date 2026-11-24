#!/bin/bash
# =============================================================================
# 01_build_signature_matrix.sh
# Build the cHL-specific CIBERSORTx signature matrix from the annotated
# scRNA-seq + snRNA-seq reference (16 cell types, 44,769 cells, 4,000 feature genes).
#
# Manuscript parameters (Methods, "Construction and benchmarking of the cHL-specific
# signature matrix"): min.expression = 0, min.genes = 400, max.genes = 800,
# hematopoietic genes = TRUE, single-cell mode.
#
# Output (in $OUTDIR): CIBERSORTx_<refsample>_inferred_phenoclasses.CIBERSORTx_<refsample>_inferred_refsample.bm.K999.txt
#   -> this is the signature matrix shipped here as resources/cHL_signature_matrix_400_800.txt
#
# Usage:  bash 01_build_signature_matrix.sh
# Requires: CIBERSORTx Fractions Singularity image + registered credentials
#           (https://cibersortx.stanford.edu) exported as CIBERSORTX_USER / CIBERSORTX_TOKEN.
# =============================================================================

CIBERSORT_ROOT="${CIBERSORT_ROOT:-/oak/stanford/groups/emoding/analysis/brian/Cibersort}"
INPUT_DIR="${INPUT_DIR:-$CIBERSORT_ROOT/bulk_input}"          # bound to /src/data inside the container
FRACTIONS_SIF="${FRACTIONS_SIF:-$CIBERSORT_ROOT/fractions.sif}"
PARTITION="${SLURM_PARTITION:-emoding}"

# Reference sample = tab-delimited count matrix, genes x cells; header row = cell-type label per cell.
# (Exported from the annotated Seurat object with the 4,000 variable features.)
refsample="/src/data/cHL_scRNA_reference_4000.txt"
outdir="$CIBERSORT_ROOT/cHL_sigmtx_400_800/"

: "${CIBERSORTX_USER:?export CIBERSORTX_USER (registered CIBERSORTx e-mail)}"
: "${CIBERSORTX_TOKEN:?export CIBERSORTX_TOKEN (token from cibersortx.stanford.edu)}"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_sigmtx
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=${PARTITION}

mkdir -p ${outdir}
singularity exec -c -B ${INPUT_DIR}:/src/data -B ${outdir}:/src/outdir \\
  ${FRACTIONS_SIF} /src/CIBERSORTxFractions \\
  --username ${CIBERSORTX_USER} --token ${CIBERSORTX_TOKEN} \\
  --refsample ${refsample} \\
  --single_cell TRUE \\
  --fraction 0 \\
  --filter TRUE \\
  --G.min 400 \\
  --G.max 800
EOF
