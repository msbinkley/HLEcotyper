# =============================================================================
# v2_G6_exvivo_recovery_inputs.R — inputs for the EcoTyper scRNA recovery of the FINAL round-14 model in the ex vivo
# Flex snRNA samples (Fig 8e); author 2026-10-25.
#
# Same as v1 (Ecotyper/ecotyper-master/2026_H2/2026-10-12_ecotyper_recovery_exvivo_inputs.R), except the labels:
#   - "B" cells take their label from the B sub-clustering (v2_G5_exvivo_B_subcluster.R, 2026-10-25; author):
#       cluster 0 -> B_postGC (class-switched memory), cluster 1 -> B_preGC (naive);
#       clusters 2 (myeloid-contaminated, 76 cells) and 3 (stress, 32 cells) are left out
#   - NK -> "NK" (the final model's name; v1 "NK_cell")
#   - all other types as v1; Innate_T, IEG_stress and Endothelial left out (as v1)
# Files (Ecotyper/datasets/):
#   2026-10-25_cHL_r14_exvivo_count.txt        genes x cells, raw UMI, first column "Gene"; cell IDs "-" -> "."
#   2026-10-25_cHL_r14_exvivo_annotation.txt   ID, CellType, Sample, Condition
#
# Run on Sherlock (in the job of v2_G6_recovery_exvivo.sh):  Rscript v2_G6_exvivo_recovery_inputs.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# Step 1. Files
F6_FILE      <- "${PROJECT_ROOT}/scRNA/2026-10-04_cHL_ex_vivo/output/2026-10-04_chl_exvivo_final_seurat_obj_f6.rds"
B_CELLS_FILE <- "${PROJECT_ROOT}/cHL_v2/output/stepE/exvivo_B/2026-10-25_exvivo_B_cells.tsv.gz"
COUNT_FILE   <- "${PROJECT_ROOT}/Ecotyper/datasets/2026-10-25_cHL_r14_exvivo_count.txt"
ANNOT_FILE   <- "${PROJECT_ROOT}/Ecotyper/datasets/2026-10-25_cHL_r14_exvivo_annotation.txt"
GENE_BLOCK   <- 1000                     # genes written per block (memory cap)

# Step 2. Ex vivo object and the B sub-cluster of each B cell
exvivo <- readRDS(F6_FILE)
b_cells <- fread(cmd = sprintf("gzip -dc %s", B_CELLS_FILE))
b_cluster_to_label <- c("0" = "B_postGC", "1" = "B_preGC")          # clusters 2 and 3 left out (author 2026-10-25)

# Step 3. Our cell types -> the final model's names (v1 map, with NK -> NK; B handled in step 4)
type_to_model <- c(
  CD4_T      = "T_CD4",
  Treg       = "T_Treg",
  CD8_T      = "T_CD8",
  NK         = "NK",
  Plasma     = "Plasma",
  Macrophage = "Macrophage",
  Monocyte   = "Macrophage",
  cDC        = "DC_Conventional",
  Fibroblast = "CAF",
  HRS        = "HRS")
model_type <- unname(type_to_model[exvivo$cell_type])                 # NA for B and the left-out types

# Step 4. B cells: label from the sub-cluster (NA for clusters 2 and 3)
is_b <- exvivo$cell_type == "B"
stopifnot(setequal(colnames(exvivo)[is_b], b_cells$cell))
b_label <- unname(b_cluster_to_label[as.character(b_cells$cluster[match(colnames(exvivo)[is_b], b_cells$cell)])])
model_type[is_b] <- b_label
keep <- !is.na(model_type)
cat(sprintf("cells kept: %d of %d (B: %d of %d)\n", sum(keep), length(keep), sum(keep & is_b), sum(is_b)))

# Step 5. Annotation (cell IDs with "-" -> ".", as EcoTyper reads them)
cell_ids <- gsub("-", ".", colnames(exvivo)[keep])
annotation <- data.table(ID = cell_ids, CellType = model_type[keep], Sample = exvivo$sample[keep],
                         Condition = as.character(exvivo$condition[keep]))
print(table(annotation$CellType, annotation$Sample))
fwrite(annotation, ANNOT_FILE, sep = "\t")

# Step 6. Counts: header line, then blocks of GENE_BLOCK genes (dense only within a block)
counts <- LayerData(exvivo, assay = "RNA", layer = "counts")[, keep]
rm(exvivo)
invisible(gc())
writeLines(paste(c("Gene", cell_ids), collapse = "\t"), COUNT_FILE)
for (block_start in seq(1, nrow(counts), by = GENE_BLOCK)) {
  block_rows <- block_start:min(block_start + GENE_BLOCK - 1, nrow(counts))
  block <- as.matrix(counts[block_rows, , drop = FALSE])
  fwrite(data.table(Gene = rownames(block), block), COUNT_FILE, sep = "\t", append = TRUE, col.names = FALSE)
}
cat(sprintf("wrote %s (%d genes x %d cells) and %s\n", COUNT_FILE, nrow(counts), ncol(counts), ANNOT_FILE))
