#!/bin/bash
# =============================================================================
# run_cytospace.sh — map single cells onto Visium spots for the four cHL tumors.
#
# Inputs (see README for how each is generated):
#   scRNA_data.txt        genes x cells raw counts from the annotated Seurat object
#                         (generate_cytospace_from_seurat_object.R)
#   cell_type_labels.txt  "Cell IDs" \t "CellType"   (same script)
#   ST_data_<tumor>.csv   genes x spots counts from Space Ranger
#   Coordinates_<tumor>.csv  spot row/col        (generate_cytospace_input_from_spaceranger_output.R)
#
# Outputs (per tumor, in $OUT/<tumor>_results): assigned_locations.csv,
#   cell_type_assignments_by_spot.csv, fractional_abundances_by_spot.csv, plots.
#   fractional_abundances_by_spot.csv (renormalised to 1 per spot) is the "Recovery cell type
#   fractions" for the *_Cyto EcoTyper Visium recovery (02_ecotyper/06_recovery_visium.sh).
#
# Usage:  bash run_cytospace.sh
# =============================================================================
CYTOSPACE_ROOT="${CYTOSPACE_ROOT:-/path/to/Cytospace}"
SC="$CYTOSPACE_ROOT/cHL_scRNA_Input"
ST="$CYTOSPACE_ROOT/cHL_Visium_Input"
OUT="$CYTOSPACE_ROOT/output"
PARTITION="${SLURM_PARTITION:-emoding}"

for t in 18081 42988 42989 47166; do
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=Cyto_${t}
#SBATCH --time=0-02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=${PARTITION}

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cytospace
mkdir -p ${OUT}/${t}_results
cd ${OUT}/${t}_results
cytospace \\
   --scRNA-path        ${SC}/scRNA_data.txt \\
   --cell-type-path    ${SC}/cell_type_labels.txt \\
   --st-path           ${ST}/ST_data_${t}.csv \\
   --coordinates-path  ${ST}/Coordinates_${t}.csv \\
   --solver-method     lap_CSPR
EOF
done
