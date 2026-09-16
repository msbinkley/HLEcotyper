#!/bin/bash
#SBATCH --job-name=Ecotyper_recovery_scRNA
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=200G
#SBATCH --partition=emoding
# =============================================================================
# 05_recovery_scRNA_external.sh — recover cell states/ecotypes in single-cell data.
#
# (a) Three external cHL scRNA-seq datasets (34 cases) concatenated into one count matrix
#     ("mega"), with an annotation file whose `Dataset` column is used for -c.
#     -t threads; -z FALSE (no z-scoring); -s 10000 = subsample cells per cell type.
# (b) Our own scRNA-seq of the four Visium tumors, one run per tumor (-c Sample), used to
#     compare in-situ (Visium) and dissociated recovery of the same tumor.
#
# Count matrices: genes x cells (raw UMI counts). Annotation: ID, CellType, <-c column>.
# Cell-type labels must use the 13 pooled names (see 01_cibersortx/resources/cell_type_mapping_16_to_13.tsv).
# =============================================================================
module load R/4.2.0
DATASETS="${DATASETS:-/path/to/datasets}"

# (a) external datasets
Rscript EcoTyper_recovery_scRNA.R -d cHL_discovery \
    -m ${DATASETS}/external_scRNA_count.txt \
    -a ${DATASETS}/external_scRNA_annotation.txt \
    -o external_scRNA_recovery \
    -c Dataset \
    -t 14 \
    -z FALSE \
    -s 10000

# (b) in-house scRNA-seq of the Visium tumors
for t in 18081 42988 42989 47166; do
  Rscript EcoTyper_recovery_scRNA.R -d cHL_discovery \
    -m ${DATASETS}/cHL_scRNA_count_${t}.txt \
    -a ${DATASETS}/cHL_scRNA_annotation_${t}.txt \
    -o cHL_scRNA_recovery_${t} \
    -c Sample
done
