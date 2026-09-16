#!/bin/bash
# =============================================================================
# 03_benchmark_HRS_subsample.sh
# Sensitivity of signature-matrix construction to the number of HRS cells in the
# reference. The reference is re-exported with HRS cells sub-sampled to 0.5 %, 1 %, 2 %,
# 5 % ... 100 % of the original, and a signature matrix is rebuilt from each
# (same parameters as 01_build_signature_matrix.sh). Each rebuilt matrix is then
# benchmarked with 02_benchmark_pseudobulk.sh.
#
# Usage:  bash 03_benchmark_HRS_subsample.sh
# =============================================================================

CIBERSORT_ROOT="${CIBERSORT_ROOT:-/oak/stanford/groups/emoding/analysis/brian/Cibersort}"
INPUT_DIR="${INPUT_DIR:-$CIBERSORT_ROOT/bulk_input}"
FRACTIONS_SIF="${FRACTIONS_SIF:-$CIBERSORT_ROOT/fractions.sif}"
PARTITION="${SLURM_PARTITION:-normal}"
: "${CIBERSORTX_USER:?}"; : "${CIBERSORTX_TOKEN:?}"

pcts=(0.5 1 2 5 10 20 30 40 50 60 70 80 90 100)

for p in "${pcts[@]}"; do
  refsample="/src/data/SubSample_HRS_${p}pct_rep1.txt"       # reference with HRS down-sampled to p %
  outdir="$CIBERSORT_ROOT/cHL_sigmtx_HRS_${p}pct/"
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Ciber_HRS_${p}
#SBATCH --time=0-00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
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
done
