# =============================================================================
# v2_C9_pseudobulk_merged.R
#
# PURPOSE : Extra pseudo-bulk benchmark set "Merged" (author 2026-10-23): ALL round-13 cells of every dataset (lab sc/sn, S1,
#           Aoki, Stewart, 10x), because the signature matrix was built from the merged reference. Same labels and design
#           as v2_C4 (final labels; 10x generic B_cell kept as "B_cell"; v1 make_pseudobulk(): 0-40 % spike-in in 2 % steps,
#           10 bootstraps, CPM).
# INPUT   : the saved step-C counts + labels (output/stepC/cHL_v2_r10_counts_universe.rds)
# OUTPUT  : output/stepC/pseudobulk/pseudobulk_Merged/pseudobulk_Merged_<type>_cpm.txt
# RUN     : on Sherlock, submitted by v2_C9_pseudobulk_merged.sh
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

# Step 2. Pseudo-bulk label = final label, except the 10x cells published as generic "B cell" (flagged in ref_exclude),
#         which are kept as their own type "B_cell" (v1 S2D convention)
cell_table[, pb_label := fifelse(dataset == "tenx" & !is.na(ref_exclude) & grepl("^10x generic B_cell", ref_exclude),
                                 "B_cell", label)]
set_cells <- cell_table[!is.na(pb_label)]                         # every labelled cell of every dataset

# Step 3. One mixture file per cell type (spike-in of that type into all other cells)
set_dir <- file.path(OUT, "pseudobulk", "pseudobulk_Merged")
dir.create(set_dir, showWarnings = FALSE, recursive = TRUE)
set_counts <- cnt[, set_cells$cell]
type_sizes <- table(set_cells$pb_label)
for (cell_type in names(type_sizes)) {
  mixtures <- make_pseudobulk(set_counts, set_cells$pb_label, cell_type)
  write_pseudobulk(mixtures, file.path(set_dir, sprintf("pseudobulk_Merged_%s_cpm.txt", cell_type)))
}
cat(sprintf("pseudo-bulk Merged: %d cells, %d types (%s)\n", nrow(set_cells), uniqueN(set_cells$pb_label),
            paste(names(type_sizes), type_sizes, sep = "=", collapse = " ")))
