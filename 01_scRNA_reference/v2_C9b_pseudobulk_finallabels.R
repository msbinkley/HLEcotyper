# =============================================================================
# v2_C9b_pseudobulk_finallabels.R
#
# PURPOSE : Author 2026-10-23: in the merged round-13 object the 157 10x cells published as generic "B cell" carry final
#           (cluster) labels (B_preGC 121, B_Memory 35, T_CD4 1), so the pseudo-bulk sets use the FINAL label for them too
#           (no "B_cell" type; v2_C4 / v2_C9 had kept v1's S2D convention). Rebuilds the two affected sets, 10X and Merged,
#           with v1 make_pseudobulk() (0-40 %, 2 % steps, 10 bootstraps, CPM).
#           Old files: output/stepC/pseudobulk/superseded_B_cell_label/ (moved there by v2_C9b_pseudobulk_finallabels.sh).
# INPUT   : output/stepC/cHL_v2_r10_counts_universe.rds (step-C counts + labels)
# OUTPUT  : output/stepC/pseudobulk/pseudobulk_{10X,Merged}/pseudobulk_<set>_<type>_cpm.txt
# RUN     : on Sherlock, submitted by v2_C9b_pseudobulk_finallabels.sh
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
W   <- "${PROJECT_ROOT}/cHL_v2"      # project root on OAK
OUT <- file.path(W, "output/stepC")
source(file.path(W, "script/utils_scRNA.R"))                     # v1 helpers: make_pseudobulk(), write_pseudobulk()

# Step 1. Read the step-C counts (genes x cells) and the cell table
saved      <- readRDS(file.path(OUT, "cHL_v2_r10_counts_universe.rds"))
cnt        <- saved$counts
cell_table <- data.table(saved$meta)
rm(saved)
invisible(gc())

# Step 2. Pseudo-bulk label = final label for every cell, incl. the 10x generic B cells
cell_table[, pb_label := label]

# Step 3. The two sets to rebuild: 10X alone, and Merged (all datasets)
SETS <- list(`10X` = "tenx", Merged = c("lab", "s1", "aoki", "tenx", "stewart"))
for (set_name in names(SETS)) {
  set_cells <- cell_table[dataset %in% SETS[[set_name]] & !is.na(pb_label)]
  set_dir <- file.path(OUT, "pseudobulk", paste0("pseudobulk_", set_name))
  dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)
  set_counts <- cnt[, set_cells$cell]
  type_sizes <- table(set_cells$pb_label)
  # one mixture file per cell type (spike-in of that type into all other cells of the set)
  for (cell_type in names(type_sizes)) {
    mixtures <- make_pseudobulk(set_counts, set_cells$pb_label, cell_type)
    write_pseudobulk(mixtures, file.path(set_dir, sprintf("pseudobulk_%s_%s_cpm.txt", set_name, cell_type)))
  }
  cat(sprintf("pseudo-bulk %s: %d cells, %d types (%s)\n", set_name, nrow(set_cells), uniqueN(set_cells$pb_label),
              paste(names(type_sizes), type_sizes, sep = "=", collapse = " ")))
  rm(set_counts)
  invisible(gc())
}
