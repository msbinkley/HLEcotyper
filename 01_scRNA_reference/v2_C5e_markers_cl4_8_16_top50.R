# =============================================================================
# v2_C5e_markers_cl4_8_16_top50.R
#
# Purpose: marker genes for the r11_rolling cluster 4 / 8 / 16 sub-UMAP (made by v2_C5c) at one resolution,
#   with the author's convention (NLPHL find_markers):
#   FindAllMarkers on log-normalised RNA, only.pos, min.pct 0.25, logfc.threshold 0.25, no downsampling;
#   top 50 per subcluster ranked by p_val_adj, then avg_log2FC; report how many of the top 50 have p_val_adj < 0.05 (QC: 50/50).
#   (The first marker table, in v2_C5c, used min.pct 0.1, 500 cells per subcluster and a log2FC ranking,
#    which favoured rare genes.)
#   Also writes cell -> subcluster for that resolution (for the per-subcluster QC table).
#
# Usage:  Rscript v2_C5e_markers_cl4_8_16_top50.R <sub-UMAP obj.rds> <out_dir> [resolution, default 0.8]
# Output (in <out_dir>):
#   2026-10-15_r11rolling_sub_cl4_8_16_cells_res<res>.tsv.gz
#   2026-10-15_r11rolling_sub_cl4_8_16_top50_markers_res<res>.tsv
#   2026-10-15_r11rolling_sub_cl4_8_16_top50_markers_res<res>_summary.tsv
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})
args    <- commandArgs(TRUE)
sub     <- readRDS(args[1])
out_dir <- args[2]
res     <- if (length(args) >= 3) args[3] else "0.8"      # clustering resolution (default 0.8)


# Step 1. Use the subclusters at this resolution; save cell -> subcluster
DefaultAssay(sub) <- "RNA"
Idents(sub) <- paste0("RNA_snn_res.", res)
cell_subcluster <- data.table(cell = colnames(sub), subcluster = as.character(Idents(sub)))
fwrite(cell_subcluster, file.path(out_dir, sprintf("2026-10-15_r11rolling_sub_cl4_8_16_cells_res%s.tsv.gz", res)), sep = "\t")


# Step 2. Markers (author rule): top 50 per subcluster, ranked by p_val_adj, then by avg_log2FC (largest first)
markers <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
top50 <- markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50)
fwrite(top50, file.path(out_dir, sprintf("2026-10-15_r11rolling_sub_cl4_8_16_top50_markers_res%s.tsv", res)), sep = "\t")


# Step 3. QC summary per subcluster: genes in the top 50, how many have p_adj < 0.05, and the first 15 genes
subclusters  <- unique(top50$cluster)                                 # in subcluster order
summary_rows <- list()
for (i in seq_along(subclusters)) {
  subcluster_markers <- top50[top50$cluster == subclusters[i], ]
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    cluster = subclusters[i],
    n_top50 = nrow(subcluster_markers),
    n_sig   = sum(subcluster_markers$p_val_adj < 0.05),
    top15   = paste(head(subcluster_markers$gene, 15), collapse = " ")
  )
}
summary_table <- rbindlist(summary_rows)
fwrite(summary_table, file.path(out_dir, sprintf("2026-10-15_r11rolling_sub_cl4_8_16_top50_markers_res%s_summary.tsv", res)), sep = "\t")
print(summary_table[order(as.integer(as.character(cluster)))])
