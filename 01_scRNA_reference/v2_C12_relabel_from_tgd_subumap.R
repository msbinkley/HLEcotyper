# =============================================================================
# v2_C12_relabel_from_tgd_subumap.R
#
# Purpose (author, 2026-10-23):
#   Re-annotate the T_CD8 / NK / ILC / T_gd cells of the merged cHL sc/sn reference from the sub-UMAP made by
#   v2_C11_subumap_r13_cd8_nk_ilc_tgd.R (round-13 clusters 3, 8, 14; 13,026 cells), at resolution 0.4:
#     subcluster 7      -> T_gd
#     subclusters 2, 10 -> NK
#     subcluster 6      -> ILC
#     all other subclusters (0, 1, 3, 4, 5, 8, 9, 11) -> T_CD8
#   Every cell in the sub-UMAP takes the label of its subcluster (the old per-cell T_gd label is no longer used).
#   Cells outside the sub-UMAP keep their round-13 label, except the 11 T_gd cells that sit in the T_CD4 cluster:
#   they are handled by OUTSIDE_TGD (see Step 3).
#
# Inputs (Box): data_v2/02_intermediate/cpA/2026-10-21_round13_labels.tsv.gz
#               data_v2/02_intermediate/subUMAP/r13_cd8_nk_ilc_tgd/2026-10-23_r13_sub_cd8_nk_ilc_tgd_cells.tsv.gz
# Outputs (Box): data_v2/02_intermediate/cpA/2026-10-23_round13_labels_tgd_subumap.tsv.gz  (same columns as the round-13 file;
#                 final_label updated; rows of removed cells dropped)
#                data_v2/02_intermediate/cpA/2026-10-23_round13_labels_tgd_subumap_changes.tsv (old label x new label counts)
# Run on the Mac:  Rscript v2_C12_relabel_from_tgd_subumap.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))
v2_dir <- "${BOX_ROOT}/paper/Manuscript/data_v2/02_intermediate"

# What to do with the T_gd cells that are not in the sub-UMAP (author decision):
#   "remove" = drop them from the reference; "T_CD4" = give them their cluster's type
OUTSIDE_TGD <- Sys.getenv("OUTSIDE_TGD", "remove")
stopifnot(OUTSIDE_TGD %in% c("remove", "T_CD4"))


# Step 1. Read the round-13 labels and the sub-UMAP cells
labels <- fread(cmd = sprintf("gzip -dc '%s'", file.path(v2_dir, "cpA/2026-10-21_round13_labels.tsv.gz")))
sub_cells <- fread(cmd = sprintf("gzip -dc '%s'",
                   file.path(v2_dir, "subUMAP/r13_cd8_nk_ilc_tgd/2026-10-23_r13_sub_cd8_nk_ilc_tgd_cells.tsv.gz")))
stopifnot(all(sub_cells$cell %in% labels$cell))
cat("round-13 cells:", nrow(labels), "; sub-UMAP cells:", nrow(sub_cells), "\n")


# Step 2. New label for every sub-UMAP cell, from its res-0.4 subcluster
subcluster_to_label <- c("7" = "T_gd",
                         "2" = "NK", "10" = "NK",
                         "6" = "ILC",
                         "0" = "T_CD8", "1" = "T_CD8", "3" = "T_CD8", "4" = "T_CD8",
                         "5" = "T_CD8", "8" = "T_CD8", "9" = "T_CD8", "11" = "T_CD8")
stopifnot(all(sub_cells$sub_res0.4 %in% names(subcluster_to_label)))
sub_cells$new_label <- subcluster_to_label[as.character(sub_cells$sub_res0.4)]

labels$old_label <- labels$final_label
row_in_labels <- match(sub_cells$cell, labels$cell)
labels$final_label[row_in_labels] <- sub_cells$new_label


# Step 3. T_gd cells outside the sub-UMAP (all in the T_CD4 cluster)
outside_tgd <- labels$old_label == "T_gd" & !(labels$cell %in% sub_cells$cell)
cat("T_gd cells outside the sub-UMAP:", sum(outside_tgd), "-> ", OUTSIDE_TGD, "\n")
if (OUTSIDE_TGD == "T_CD4") {
  labels$final_label[outside_tgd] <- "T_CD4"
} else {
  labels$final_label[outside_tgd] <- "REMOVED"
}


# Step 4. Table of changes (old label x new label), for the log
changes <- labels[old_label != final_label, .N, by = .(old_label, new_label = final_label)]
changes <- changes[order(-N)]
print(changes)
print(labels[, .(before = sum(old_label == "T_gd"), after = sum(final_label == "T_gd"))])
fwrite(changes, file.path(v2_dir, "cpA/2026-10-23_round13_labels_tgd_subumap_changes.tsv"), sep = "\t")


# Step 5. Write the new label file (removed cells dropped; same columns as the round-13 file)
kept <- labels[final_label != "REMOVED"]
kept$old_label <- NULL
out_file <- file.path(v2_dir, "cpA/2026-10-23_round13_labels_tgd_subumap.tsv")
fwrite(kept, out_file, sep = "\t")
system(sprintf("gzip -f '%s'", out_file))
cat("written:", paste0(out_file, ".gz"), "-", nrow(kept), "cells\n")
print(kept[, .N, by = final_label][order(-N)])
