# =============================================================================
# v2_C5h_sub3_tgd_vs_nk.R
#
# Purpose (author 2026-10-16): split round-4 sub-UMAP subcluster 3 (NK CD16+, 736 cells, 278 published T_gd) into T_gd and NK.
#   Author: "look at T_gd scores using the 278 published T_gd cells in 3 as reference; cells close to 14 seem to have
#   higher T_gd genes; can we run a sub-UMAP of 3 and break it into T_gd and NK?"
#
# Steps:
#   1. T_gd-vs-NK signature from published labels inside subcluster 3 (non-S1 cells): FindMarkers T_gd vs NK, both directions,
#      author marker rule (p_adj-ranked, p_adj < 0.05, min.pct 0.25, logfc 0.25), top 30 genes each way.
#      score = AddModuleScore(T_gd genes) - AddModuleScore(NK genes); AUC of the score for published T_gd vs NK.
#   2. Score plotted on the round-4 sub-UMAP (all cells) to check that cells near subcluster 14 score higher.
#   3. Sub-UMAP of subcluster 3 alone: LogNormalize -> 2,000 HVG (all cells; per-dataset HVG would be unstable with
#      13 Stewart cells) -> PCA 20 -> Harmony (dataset + chemistry + sample, theta 4/4/2) -> Louvain res 0.3 / 0.5 / 0.8 -> UMAP.
#      Per sub-subcluster: published T_gd / NK %, S1 cells, median score; top-50 markers (author rule) at res 0.5.
#
# Usage: Rscript v2_C5h_... <round-4 obj.rds> <Box table dir> <Box figure dir>   (runs locally)
# Output: <Box table dir>/2026-10-16_r11rolling_sub_cl4_8_16_r4_sub3_{T_gd_vs_NK_markers.tsv, composition.tsv,
#                                                                     top50_markers_res0.5.tsv, cells.tsv.gz}
#         <Box figure dir>/2026-10-16_r11rolling_sub_cl4_8_16_r4_sub3_tgd_vs_nk.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
args       <- commandArgs(TRUE)
round4     <- readRDS(args[1])       # round-4 sub-UMAP object (v2_C5f, TAG=r4)
table_dir  <- args[2]
figure_dir <- args[3]
stamp <- "2026-10-16"
tag   <- "r11rolling_sub_cl4_8_16_r4_sub3"
table_file  <- function(name) file.path(table_dir,  paste0(stamp, "_", tag, "_", name))
figure_file <- function(name) file.path(figure_dir, paste0(stamp, "_", tag, "_", name))

# Published label; S1 cells (no published label) are called "S1"
round4$pub <- ifelse(round4$dataset == "s1", "S1", round4$label)
cells_in_3 <- colnames(round4)[round4$RNA_snn_res.1.2 == "3"]          # round-4 res-1.2 subcluster 3
cat("subcluster 3:", length(cells_in_3), "cells\n")
print(table(round4$pub[cells_in_3]))


# Step 1. Signature from published T_gd vs NK inside subcluster 3
reference_cells <- cells_in_3[round4$pub[cells_in_3] %in% c("T_gd", "NK")]
reference <- subset(round4, cells = reference_cells)
Idents(reference) <- "pub"
up_in_tgd <- FindMarkers(reference, ident.1 = "T_gd", ident.2 = "NK", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
up_in_nk  <- FindMarkers(reference, ident.1 = "NK", ident.2 = "T_gd", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

# Helper: genes with p_adj < 0.05, ranked by p_adj then avg_log2FC (largest first); the first 30
top30_genes <- function(marker_table) {
  marker_table <- marker_table[marker_table$p_val_adj < 0.05, ]
  marker_table <- marker_table[order(marker_table$p_val_adj, -marker_table$avg_log2FC), ]
  head(rownames(marker_table), 30)
}
tgd_genes <- top30_genes(up_in_tgd)
nk_genes  <- top30_genes(up_in_nk)
cat("T_gd genes (", length(tgd_genes), "):", paste(tgd_genes, collapse = " "),
    "\nNK genes (", length(nk_genes), "):", paste(nk_genes, collapse = " "), "\n")

# Save both marker lists (p_adj < 0.05), sorted by direction, p_adj, then avg_log2FC (largest first)
signature <- rbind(data.table(direction = "up_in_T_gd", gene = rownames(up_in_tgd), up_in_tgd),
                   data.table(direction = "up_in_NK",   gene = rownames(up_in_nk),  up_in_nk))
signature <- signature[p_val_adj < 0.05]
signature <- signature[order(direction, p_val_adj, -avg_log2FC)]
fwrite(signature, table_file("T_gd_vs_NK_markers.tsv"), sep = "\t")

# Score = module score of the T_gd genes - module score of the NK genes (seed 42)
set.seed(42)
round4 <- AddModuleScore(round4, features = list(tgd_genes, nk_genes), name = "sigTN", seed = 42)   # adds sigTN1 (T_gd), sigTN2 (NK)
round4$Tgd_score <- round4$sigTN1 - round4$sigTN2

# Helper: AUC = probability that a random case scores higher than a random control (Mann-Whitney, from ranks)
auc <- function(case_scores, control_scores) {
  ranks <- rank(c(case_scores, control_scores))
  n_case    <- length(case_scores)
  n_control <- length(control_scores)
  (sum(ranks[seq_along(case_scores)]) - n_case * (n_case + 1) / 2) / (n_case * n_control)
}
score_in_3 <- round4$Tgd_score[cells_in_3]
pub_in_3   <- round4$pub[cells_in_3]
cat(sprintf("AUC (published T_gd vs NK, in subcluster 3): %.3f\n", auc(score_in_3[pub_in_3 == "T_gd"], score_in_3[pub_in_3 == "NK"])))
score_table <- data.table(pub = pub_in_3, score = score_in_3)
print(score_table[, .(n = .N, median_score = round(median(score), 3)), by = pub])


# Step 2. Score on the round-4 sub-UMAP (all cells, then subclusters 3 and 14 only)
umap_table <- data.table(Embeddings(round4, "umap"), keep.rownames = "cell")
umap_table[, score := round4$Tgd_score]
umap_table[, sub := as.character(round4$RNA_snn_res.1.2)]
umap_table[, pub := round4$pub]
setnames(umap_table, 2:3, c("U1", "U2"))

all_by_score   <- umap_table[order(score)]                           # highest scores drawn on top
in_3_14        <- umap_table[sub %in% c("3", "14")]
in_3_14_by_score <- in_3_14[order(score)]

score_all_plot <- ggplot(all_by_score, aes(U1, U2, colour = score)) +
                  geom_point(size = 0.6, stroke = 0) +
                  coord_equal() +
                  theme_classic(12) +
                  scale_colour_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0, name = "T_gd - NK\nscore") +
                  ggtitle("T_gd-vs-NK score on the round-4 sub-UMAP")
score_3_14_plot <- ggplot(in_3_14_by_score, aes(U1, U2, colour = score)) +
                   geom_point(size = 1.2, stroke = 0) +
                   coord_equal() +
                   theme_classic(12) +
                   scale_colour_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0, name = "T_gd - NK\nscore") +
                   ggtitle("subclusters 3 and 14 only")
published_3_14_plot <- ggplot(in_3_14, aes(U1, U2, colour = pub)) +
                       geom_point(size = 1.2, stroke = 0, alpha = 0.8) +
                       coord_equal() +
                       theme_classic(12) +
                       scale_colour_manual(values = c(T_gd = "#A65628", NK = "#FDBF6F", T_CD8 = "#999999", S1 = "#5AAE61", ILC = "#00C19A"),
                                           name = "published") +
                       ggtitle("subclusters 3 and 14: published labels")


# Step 3. Sub-UMAP of subcluster 3 alone
#   LogNormalize -> 2,000 HVG over all cells -> scale -> PCA 20 -> Harmony (theta 4/4/2) -> Louvain res 0.3 / 0.5 / 0.8
#   -> UMAP (20 neighbours); seeds 42
sub3_counts <- LayerData(round4, assay = "RNA", layer = "counts")[, cells_in_3]
sub3_meta   <- round4@meta.data[cells_in_3, c("dataset", "sample", "chemistry", "label", "pub", "Tgd_score")]
sub3 <- CreateSeuratObject(sub3_counts, meta.data = sub3_meta)
sub3 <- NormalizeData(sub3, verbose = FALSE)
sub3 <- FindVariableFeatures(sub3, nfeatures = 2000, verbose = FALSE)
sub3 <- ScaleData(sub3, verbose = FALSE)                              # scales the 2,000 variable genes
sub3 <- RunPCA(sub3, npcs = 20, seed.use = 42, verbose = FALSE)
set.seed(42)
sub3 <- RunHarmony(sub3, group.by.vars = c("dataset", "chemistry", "sample"), theta = c(4, 4, 2),
                   reduction.use = "pca", reduction.save = "harmony", verbose = FALSE)
sub3 <- FindNeighbors(sub3, reduction = "harmony", dims = 1:20, verbose = FALSE)
for (res in c(0.3, 0.5, 0.8)) {
  sub3 <- FindClusters(sub3, resolution = res, random.seed = 42, verbose = FALSE)   # adds RNA_snn_res.<res>
}
sub3 <- RunUMAP(sub3, reduction = "harmony", dims = 1:20, seed.use = 42, n.neighbors = 20, verbose = FALSE)

# Composition of every sub-subcluster: published T_gd / NK / S1 / other cells, % T_gd among T_gd + NK, median score
composition <- list()
for (res_value in c("0.3", "0.5", "0.8")) {
  cell_info <- data.table(cl    = as.character(sub3@meta.data[[paste0("RNA_snn_res.", res_value)]]),
                          pub   = sub3$pub,
                          score = sub3$Tgd_score)
  res_table <- cell_info[, .(
      cells  = .N,
      T_gd   = sum(pub == "T_gd"),
      NK     = sum(pub == "NK"),
      S1     = sum(pub == "S1"),
      other  = sum(!pub %in% c("T_gd", "NK", "S1")),
      pct_T_gd_of_labelled = round(100 * sum(pub == "T_gd") / max(1, sum(pub %in% c("T_gd", "NK"))), 1),   # max(1, ...) avoids 0 / 0
      median_score = round(median(score), 3)
    ), by = cl]                                                       # sub-subclusters in order of first appearance
  res_table[, res := res_value]
  res_table <- res_table[order(-pct_T_gd_of_labelled)]
  composition[[length(composition) + 1]] <- res_table
}
composition <- rbindlist(composition)
setcolorder(composition, "res")                                       # res first
fwrite(composition, table_file("composition.tsv"), sep = "\t")
print(composition)

# Markers at res 0.5 (author rule): top 50 per sub-subcluster, ranked by p_val_adj, then by avg_log2FC (largest first)
Idents(sub3) <- "RNA_snn_res.0.5"
markers <- FindAllMarkers(sub3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
top50 <- markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50)
fwrite(top50, table_file("top50_markers_res0.5.tsv"), sep = "\t")

# Print a QC summary: how many of the top 50 have p_adj < 0.05, and the first 15 genes
subclusters  <- unique(top50$cluster)                                 # in subcluster order
summary_rows <- list()
for (i in seq_along(subclusters)) {
  subcluster_markers <- top50[top50$cluster == subclusters[i], ]
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    cluster = subclusters[i],
    n_sig   = sum(subcluster_markers$p_val_adj < 0.05),
    top15   = paste(head(subcluster_markers$gene, 15), collapse = " ")
  )
}
print(rbindlist(summary_rows))

# Cell table: sub-subclusters at the 3 resolutions, score, UMAP coordinates
cell_out <- data.table(cell        = colnames(sub3),
                       sub3_res0.3 = sub3$RNA_snn_res.0.3,
                       sub3_res0.5 = sub3$RNA_snn_res.0.5,
                       sub3_res0.8 = sub3$RNA_snn_res.0.8,
                       Tgd_score   = sub3$Tgd_score,
                       Embeddings(sub3, "umap"))                      # adds umap_1, umap_2
fwrite(cell_out, table_file("cells.tsv.gz"), sep = "\t")


# Step 4. Figures
published_colors <- c(T_gd = "#A65628", NK = "#FDBF6F", T_CD8 = "#999999", S1 = "#5AAE61", ILC = "#00C19A")
pdf(figure_file("tgd_vs_nk.pdf"), width = 20, height = 8)
# page 1: the score on the round-4 sub-UMAP
print(score_all_plot | score_3_14_plot | published_3_14_plot)
# page 2: sub-UMAP of subcluster 3: published label | score | sub-subclusters at res 0.5
sub3_label_plot <- DimPlot(sub3, group.by = "pub", pt.size = 1.5, cols = published_colors) +
                   ggtitle(sprintf("sub-UMAP of subcluster 3 (%d cells): published label", ncol(sub3)))
sub3_score_plot <- FeaturePlot(sub3, features = "Tgd_score", pt.size = 1.5) +
                   scale_colour_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0) +
                   ggtitle("T_gd - NK score")
sub3_res05_plot <- DimPlot(sub3, group.by = "RNA_snn_res.0.5", pt.size = 1.5, label = TRUE) + ggtitle("sub-subclusters, res 0.5")
print(sub3_label_plot | sub3_score_plot | sub3_res05_plot)
# page 3: res 0.3 | res 0.8 | dataset
sub3_res03_plot   <- DimPlot(sub3, group.by = "RNA_snn_res.0.3", pt.size = 1.5, label = TRUE) + ggtitle("res 0.3")
sub3_res08_plot   <- DimPlot(sub3, group.by = "RNA_snn_res.0.8", pt.size = 1.5, label = TRUE) + ggtitle("res 0.8")
sub3_dataset_plot <- DimPlot(sub3, group.by = "dataset", pt.size = 1.5) + ggtitle("dataset")
print(sub3_res03_plot | sub3_res08_plot | sub3_dataset_plot)
# page 4: top 6 signature genes each way + T-cell / gamma-delta / NK genes
feature_genes <- intersect(c(head(tgd_genes, 6), head(nk_genes, 6), "CD3E", "CD3D", "TRDC", "TRGC1", "TRGC2", "FCGR3A"), rownames(sub3))
print(FeaturePlot(sub3, features = feature_genes, ncol = 6, pt.size = 0.8, order = TRUE))
dev.off()
cat("done\n")
