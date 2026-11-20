# =============================================================================
# v2_C19_record_subumap_removals_r15.R
#
# Purpose (author, 2026-10-25):
#   "record c7 in subclusters, res 0.4 as tgd. record cells that were removed in subumap, then remove these cells from
#    r14 to r15. record cells in c7 to tgd in r15. recluster r15"
#   This script writes the two records (removed cells, T_gd cells) and the round-15 label file.
#   The recluster is v2_C20_recluster_r15.R on Sherlock.
#
# The sub-UMAP chain (round-14 clusters 4 / 9 / 14 = T_CD8 / NK / ILC; all 717 round-14 T_gd cells are in it):
#   r14 sub-UMAP (v2_C15, 12,763 cells)
#     step 1: T_gd cells outside res-1.2 c12 removed (337); c12 cells labelled T_gd               -> r14b (v2_C17)
#   r14b sub-UMAP (12,426 cells)
#     step 2: res-1.2 c15 removed (187); T_gd cells outside res-1.2 c8 removed (24)               -> r14c (v2_C18)
#   r14c sub-UMAP (12,215 cells)
#     step 3: res-0.4 c7 = T_gd (357 cells); T_gd cells outside c7 removed
#             (author 2026-10-25, answer to "What should the 21 other T_gd be in r15?": "Remove them")
#
# Round 15 = round 14 minus every removed cell; label = round-14 label, except every r14c res-0.4 c7 cell = T_gd.
#
# Inputs (Box): data_v2/02_intermediate/cpA/2026-10-23_round14_labels.tsv.gz
#               data_v2/02_intermediate/subUMAP/r14_cd8_nk_ilc_tgd/2026-10-25_r14_sub_cd8_nk_ilc_tgd_cells.tsv.gz
#               data_v2/02_intermediate/subUMAP/r14b_cd8_nk_ilc_tgd/2026-10-25_r14b_sub_cd8_nk_ilc_tgd_cells.tsv.gz
#               data_v2/02_intermediate/subUMAP/r14c_cd8_nk_ilc_tgd/2026-10-25_r14c_sub_cd8_nk_ilc_tgd_cells.tsv.gz
# Outputs (Box, data_v2/02_intermediate/cpA/):
#   2026-10-25_r15_removed_cells.tsv.gz   one row per removed cell: cell, dataset, sample, r14_label, step, reason
#   2026-10-25_r15_tgd_cells.tsv.gz       one row per round-15 T_gd cell (r14c res-0.4 c7)
#   2026-10-25_r15_label_changes.tsv      round-14 label x round-15 label counts (removed = "REMOVED")
#   2026-10-25_round15_labels_pre.tsv.gz  cell, dataset, r14_label, final_label for the round-15 cells
#                                         (clusters are added by the recluster)
# Run on the Mac:  Rscript v2_C19_record_subumap_removals_r15.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

# Step 1. Folders and the author's choices
v2_dir  <- "${BOX_ROOT}/paper/Manuscript/data_v2/02_intermediate"
sub_dir <- file.path(v2_dir, "subUMAP")
r14_tgd_subcluster  <- "12"     # step 1: res-1.2 c12 of the r14 sub-UMAP
r14b_drop_subcluster <- "15"    # step 2: res-1.2 c15 of the r14b sub-UMAP (removed)
r14b_tgd_subcluster <- "8"      # step 2: res-1.2 c8 of the r14b sub-UMAP (T_gd kept only here)
r14c_tgd_subcluster <- "7"      # step 3: res-0.4 c7 of the r14c sub-UMAP (= T_gd in round 15)
read_gz <- function(path) fread(cmd = sprintf("gzip -dc '%s'", path))


# Step 2. Read the round-14 labels and the three sub-UMAP cell tables
r14_labels <- read_gz(file.path(v2_dir, "cpA/2026-10-23_round14_labels.tsv.gz"))
sub_r14  <- read_gz(file.path(sub_dir, "r14_cd8_nk_ilc_tgd/2026-10-25_r14_sub_cd8_nk_ilc_tgd_cells.tsv.gz"))
sub_r14b <- read_gz(file.path(sub_dir, "r14b_cd8_nk_ilc_tgd/2026-10-25_r14b_sub_cd8_nk_ilc_tgd_cells.tsv.gz"))
sub_r14c <- read_gz(file.path(sub_dir, "r14c_cd8_nk_ilc_tgd/2026-10-25_r14c_sub_cd8_nk_ilc_tgd_cells.tsv.gz"))
stopifnot(all(sub_r14$cell %in% r14_labels$cell))
stopifnot(all(sub_r14b$cell %in% sub_r14$cell))
stopifnot(all(sub_r14c$cell %in% sub_r14b$cell))
cat("round 14:", nrow(r14_labels), "cells; sub-UMAPs r14 / r14b / r14c:", nrow(sub_r14), "/", nrow(sub_r14b), "/", nrow(sub_r14c), "\n")


# Step 3. The removed cells of each step, with the reason
# 3a. step 1: round-14 T_gd cells outside r14 c12
r14_c12 <- sub_r14$cell[as.character(sub_r14$sub_res1.2) == r14_tgd_subcluster]
step1 <- sub_r14[final_label == "T_gd" & !(cell %in% r14_c12), cell]

# 3b. step 2: r14b c15, and r14b T_gd cells outside r14b c8
step2_c15 <- sub_r14b[as.character(sub_res1.2) == r14b_drop_subcluster, cell]
step2_tgd <- sub_r14b[final_label == "T_gd" & as.character(sub_res1.2) != r14b_tgd_subcluster, cell]

# 3c. step 3: r14c T_gd cells outside r14c res-0.4 c7
r14c_c7 <- sub_r14c$cell[as.character(sub_r14c$sub_res0.4) == r14c_tgd_subcluster]
step3 <- sub_r14c[final_label == "T_gd" & !(cell %in% r14c_c7), cell]

# 3d. check: steps 1 and 2 are exactly the cells the sub-UMAP runs dropped
stopifnot(setequal(c(step1), setdiff(sub_r14$cell, sub_r14b$cell)))
stopifnot(setequal(c(step2_c15, step2_tgd), setdiff(sub_r14b$cell, sub_r14c$cell)))

removed <- rbind(
  data.table(cell = step1,     step = "1 (r14 sub-UMAP)",  reason = "T_gd outside res-1.2 c12"),
  data.table(cell = step2_c15, step = "2 (r14b sub-UMAP)", reason = "res-1.2 c15 (CD8 GZMK+ memory, weak)"),
  data.table(cell = step2_tgd, step = "2 (r14b sub-UMAP)", reason = "T_gd outside res-1.2 c8"),
  data.table(cell = step3,     step = "3 (r14c sub-UMAP)", reason = "T_gd outside res-0.4 c7"))
stopifnot(!anyDuplicated(removed$cell))
removed <- merge(removed, sub_r14[, .(cell, dataset, sample, r14_label = final_label)], by = "cell")
print(removed[, .N, by = .(step, reason, r14_label)])
cat("removed in total:", nrow(removed), "\n")


# Step 4. The round-15 T_gd cells = r14c res-0.4 c7
tgd_cells <- sub_r14c[cell %in% r14c_c7, .(cell, dataset, sample, r14c_label = final_label, sub_res0.4, sub_res1.2)]
tgd_cells <- merge(tgd_cells, sub_r14[, .(cell, r14_label = final_label)], by = "cell")
print(tgd_cells[, .N, by = .(r14_label, r14c_label)])


# Step 5. Round-15 labels: round 14 minus the removed cells; c7 cells -> T_gd; all others keep the round-14 label
r15 <- r14_labels[!(cell %in% removed$cell), .(cell, dataset, r14_label = final_label)]
r15$final_label <- r15$r14_label
r15[cell %in% tgd_cells$cell, final_label := "T_gd"]
stopifnot(nrow(r15) == nrow(r14_labels) - nrow(removed))
stopifnot(sum(r15$final_label == "T_gd") == nrow(tgd_cells))

# check: inside the sub-UMAP, the round-15 label equals the r14c sub-UMAP label for every kept cell
in_r14c <- merge(r15, sub_r14c[, .(cell, r14c_label = final_label)], by = "cell")
stopifnot(all(in_r14c$final_label == in_r14c$r14c_label | in_r14c$cell %in% tgd_cells$cell))
cat("round 15:", nrow(r15), "cells; T_gd", sum(r15$final_label == "T_gd"), "\n")


# Step 6. Change table (round-14 label x round-15 label; removed cells = "REMOVED")
change_rows <- rbind(r15[r14_label != final_label, .(r14_label, r15_label = final_label)],
                     removed[, .(r14_label, r15_label = "REMOVED")])
changes <- change_rows[, .N, by = .(r14_label, r15_label)][order(-N)]
print(changes)


# Step 7. Write (never over an existing file)
out_files <- file.path(v2_dir, "cpA", c("2026-10-25_r15_removed_cells.tsv", "2026-10-25_r15_tgd_cells.tsv",
                                        "2026-10-25_r15_label_changes.tsv", "2026-10-25_round15_labels_pre.tsv"))
stopifnot(!any(file.exists(out_files)), !any(file.exists(paste0(out_files, ".gz"))))
fwrite(removed[order(step, cell)], out_files[1], sep = "\t")
fwrite(tgd_cells, out_files[2], sep = "\t")
fwrite(changes, out_files[3], sep = "\t")
fwrite(r15, out_files[4], sep = "\t")
for (gz_file in out_files[c(1, 2, 4)]) {
  system(sprintf("gzip '%s'", gz_file))
}
cat("written to", file.path(v2_dir, "cpA"), "\n")
