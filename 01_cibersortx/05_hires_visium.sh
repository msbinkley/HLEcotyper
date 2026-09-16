#!/bin/bash
# =============================================================================
# 05_hires_visium.sh
# CIBERSORTx HiRes (S-mode) on the four Visium tumors, spot-level.
# Mixture = Space Ranger counts per spot (genes x spots), one file per tumor.
# The per-spot Fractions-Adjusted output is renormalised to sum to 1 and used as
# "Recovery cell type fractions" for EcoTyper Visium recovery (02_ecotyper/06_recovery_visium.sh).
#
# Usage:  bash 05_hires_visium.sh
# =============================================================================

CIBERSORT_ROOT="${CIBERSORT_ROOT:-/oak/stanford/groups/emoding/analysis/brian/Cibersort}"
INPUT_DIR="${INPUT_DIR:-$CIBERSORT_ROOT/bulk_input}"
HIRES_SIF="${HIRES_SIF:-$CIBERSORT_ROOT/hires.sif}"
PARTITION="${SLURM_PARTITION:-emoding}"
: "${CIBERSORTX_USER:?}"; : "${CIBERSORTX_TOKEN:?}"

sigmatrix="/src/data/cHL_signature_matrix_400_800.txt"
refsample="/src/data/cHL_scRNA_reference_4000.txt"
genelist="/src/data/coding_gene_list.txt"

tumors=(18081 42988 42989 47166)

for t in "${tumors[@]}"; do
  mixture="/src/data/ST_data_${t}.txt"
  outdir="$CIBERSORT_ROOT/cHL_visium_${t}/"
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_visium_${t}
#SBATCH --time=0-02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
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
done
