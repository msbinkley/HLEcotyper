# =============================================================================
# v2_C5f_subumap_cl4_8_16_drop13_recluster.R
#
# Purpose (author 2026-10-15): cleanup rounds of the cluster 4 / 8 / 16 sub-UMAP:
#   remove a trash subcluster, recluster, re-annotate. One run = one round.
#
# Settings come from environment variables (defaults = round 2):
#   TAG       round name, used in every output file name                       (default r2)
#   DROP      subcluster(s) to remove, comma-separated; may be empty            (default 13)
#   DROP_COL  column of the previous round's cell table that holds DROP         (default subcluster)
#   N_EXPECT  number of cells that must be removed (safety check)               (default 253)
#   THETA     Harmony theta for dataset, chemistry, sample                      (default 2,2,2)
#   REASON    text written next to every removed cell                           (default: the round-2 reason)
#
# Rounds run with this script:
#   round 2 : TAG=r2  DROP=13    DROP_COL=subcluster  N_EXPECT=253
#             (res-1.2 subcluster 13 of round 1: heat-shock activated CD8, 25/50)
#   round 3 : TAG=r3  DROP=19    DROP_COL=sub_res1.2  N_EXPECT=104
#             (res-1.2 subcluster 19 of round 2: MT- genes lead, 31/50)
#   round 3B: TAG=r3B DROP=      N_EXPECT=0  THETA=4,4,2
#             (same cells as round 3; Harmony theta 4/4/2 chosen by the author 2026-10-16 from the theta test, v2_C5g;
#              the default theta 2 left S1 poorly mixed)
#   round 4 : TAG=r4  DROP=10,18 DROP_COL=sub_res1.2  N_EXPECT=588  THETA=4,4,2
#             (round-3B res-1.2 subclusters 10: generic genes, 27/50 sig; 18: CD4-expressing T cells, 86% S1,
#              inside the CD8 cluster)
#
# Recipe = v2_C5c: LogNormalize -> 3,000 variable genes chosen per dataset (ranked by number of datasets) -> scale -> PCA 30
#   -> Harmony (dataset + chemistry + sample, seed 42) -> Louvain res 0.4 / 0.8 / 1.2 -> UMAP (seed 42).
# Markers (author rule): FindAllMarkers on log-normalised RNA, only.pos, min.pct 0.25, logfc 0.25, no downsampling;
#   top 50 per subcluster ranked by p_val_adj, then avg_log2FC; report how many of the 50 have p_val_adj < 0.05. Resolution 1.2.
#
# Usage: Rscript v2_C5f_... <round-1 obj.rds> <round-1 cells_res1.2.tsv.gz> <work dir for the object> <Box table dir> <Box figure dir>
#   (runs locally; the object stays out of Box because it carries Aoki cell-level data)
# Output: <work dir>/r11rolling_sub_cl4_8_16_<TAG>_obj.rds
#         <Box table dir>/2026-10-15_r11rolling_sub_cl4_8_16_<TAG>_{removed_cells.tsv.gz, cells.tsv.gz,
#                                                                 top50_markers_res1.2.tsv, top50_markers_res1.2_summary.tsv}
#         <Box figure dir>/2026-10-15_r11rolling_sub_cl4_8_16_<TAG>_{umaps.pdf, markers.pdf}
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
args           <- commandArgs(TRUE)
previous_obj   <- args[1]      # previous round's sub-UMAP object
previous_cells <- args[2]      # previous round's cell -> subcluster table
work_dir       <- args[3]      # where the new object is saved (outside Box)
table_dir      <- args[4]
figure_dir     <- args[5]

stamp <- "2026-10-15"                                          # fixed date prefix for every round
tag   <- paste0("r11rolling_sub_cl4_8_16_", Sys.getenv("TAG", "r2"))
THETA <- as.numeric(strsplit(Sys.getenv("THETA", "2,2,2"), ",")[[1]])   # Harmony theta: dataset, chemistry, sample
stopifnot(length(THETA) == 3)
DROP     <- strsplit(Sys.getenv("DROP", "13"), ",")[[1]]
DROP_COL <- Sys.getenv("DROP_COL", "subcluster")
N_EXPECT <- as.integer(Sys.getenv("N_EXPECT", "253"))
REASON   <- Sys.getenv("REASON", "round-1 res-1.2 subcluster 13: heat-shock activated CD8, 25/50 sig markers, low UMI")
table_file  <- function(name) file.path(table_dir,  paste0(stamp, "_", tag, "_", name))
figure_file <- function(name) file.path(figure_dir, paste0(stamp, "_", tag, "_", name))


# Step 1. Find the cells to remove (the DROP subclusters of the previous round) and save them with the reason
previous <- readRDS(previous_obj)
previous_cell_table <- fread(previous_cells)
drop_column_values <- as.character(previous_cell_table[, get(DROP_COL)])
removed_cells <- previous_cell_table$cell[drop_column_values %in% DROP]
stopifnot(length(removed_cells) == N_EXPECT, all(removed_cells %in% colnames(previous)))
removed_table <- data.table(cell    = removed_cells,
                            dataset = previous$dataset[removed_cells],
                            label   = previous$label[removed_cells],
                            reason  = REASON)
fwrite(removed_table, table_file("removed_cells.tsv.gz"), sep = "\t")


# Step 2. New Seurat object from the raw counts of the kept cells
counts <- LayerData(previous, assay = "RNA", layer = "counts")
kept_cells <- setdiff(colnames(previous), removed_cells)
kept_meta  <- previous@meta.data[kept_cells, c("cell", "dataset", "sample", "chemistry", "label", "s1_cell_type", "whole_cluster", "label_show")]
sub <- CreateSeuratObject(counts[, kept_cells], meta.data = kept_meta)
rm(previous, counts)
invisible(gc())
cat(ncol(sub), "cells after removing", length(removed_cells), "\n")


# Step 3. LogNormalize
sub <- NormalizeData(sub, verbose = FALSE)


# Step 4. Variable genes: 3,000 per dataset, then rank genes by how many datasets pick them
#         (ties broken by the gene's median rank across datasets, then alphabetically); keep the top 3,000
cells_by_dataset <- split(colnames(sub), sub$dataset)                 # datasets in alphabetical order
hvg_per_dataset <- list()
for (ds in names(cells_by_dataset)) {
  ds_object <- CreateSeuratObject(LayerData(sub, assay = "RNA", layer = "counts")[, cells_by_dataset[[ds]]])
  ds_object <- FindVariableFeatures(ds_object, nfeatures = 3000, verbose = FALSE)
  hvg_per_dataset[[ds]] <- VariableFeatures(ds_object)
}
times_picked    <- table(unlist(hvg_per_dataset))                     # genes in alphabetical order; value = number of datasets
candidate_genes <- names(times_picked)
median_rank <- numeric(length(candidate_genes))
for (i in seq_along(candidate_genes)) {
  ranks <- c()
  for (ds in names(hvg_per_dataset)) {
    ranks <- c(ranks, match(candidate_genes[i], hvg_per_dataset[[ds]]))   # NA if this dataset did not pick the gene
  }
  median_rank[i] <- median(ranks, na.rm = TRUE)
}
hvg <- candidate_genes[order(-times_picked, median_rank)][1:3000]
VariableFeatures(sub) <- hvg


# Step 5. Scale -> PCA 30 -> Harmony (theta from THETA) -> neighbours -> Louvain at 3 resolutions -> UMAP (all seeds 42)
sub <- ScaleData(sub, features = hvg, verbose = FALSE)
sub <- RunPCA(sub, features = hvg, npcs = 30, seed.use = 42, verbose = FALSE)
set.seed(42)
sub <- RunHarmony(sub, group.by.vars = c("dataset", "chemistry", "sample"), theta = THETA,
                  reduction.use = "pca", reduction.save = "harmony", verbose = FALSE)
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (res in c(0.4, 0.8, 1.2)) {
  sub <- FindClusters(sub, resolution = res, random.seed = 42, verbose = FALSE)   # adds RNA_snn_res.<res>
}
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)


# Step 6. Final annotation: the whole-object (r11_rolling) cluster's type,
#         except protected types keep their own published label (T_gd, ILC, HRS; non-S1 cells only)
cluster_type      <- c(`4` = "T_CD8", `8` = "NK", `16` = "ILC")
protected_types   <- c("T_gd", "ILC", "HRS")
is_protected      <- sub$dataset != "s1" & sub$label %in% protected_types
type_from_cluster <- unname(cluster_type[sub$whole_cluster])
sub$final_label <- ifelse(is_protected, sub$label, type_from_cluster)
saveRDS(sub, file.path(work_dir, paste0(tag, "_obj.rds")))

# Cell table: subclusters at the 3 resolutions + UMAP coordinates
cell_out <- data.table(cell       = colnames(sub),
                       sub_res0.4 = sub$RNA_snn_res.0.4,
                       sub_res0.8 = sub$RNA_snn_res.0.8,
                       sub_res1.2 = sub$RNA_snn_res.1.2,
                       Embeddings(sub, "umap"))                       # adds umap_1, umap_2
fwrite(cell_out, table_file("cells.tsv.gz"), sep = "\t")


# Step 7. Markers at res 1.2 (author rule): top 50 per subcluster, ranked by p_val_adj, then by avg_log2FC (largest first)
Idents(sub) <- "RNA_snn_res.1.2"
markers <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
top50 <- markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50)
fwrite(top50, table_file("top50_markers_res1.2.tsv"), sep = "\t")

# QC summary per subcluster: how many of the top 50 have p_adj < 0.05, the largest p_adj among them, the first 15 genes
subclusters  <- unique(top50$cluster)                                 # in subcluster order
summary_rows <- list()
for (i in seq_along(subclusters)) {
  subcluster_markers <- top50[top50$cluster == subclusters[i], ]
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    cluster  = subclusters[i],
    n_sig    = sum(subcluster_markers$p_val_adj < 0.05),
    max_padj = signif(max(subcluster_markers$p_val_adj), 2),
    top15    = paste(head(subcluster_markers$gene, 15), collapse = " ")
  )
}
summary_table <- rbindlist(summary_rows)
fwrite(summary_table, table_file("top50_markers_res1.2_summary.tsv"), sep = "\t")
print(summary_table[order(as.integer(as.character(cluster)))])


# Step 8. Figures
label_colors <- c(ILC = "#00C19A", NK = "#FDBF6F", T_CD8 = "#999999", T_gd = "#A65628", T_CD4 = "#B2DF8A",
                  T_Treg = "#A6CEE3", B_preGC = "#1F78B4", `S1 (unlabelled here)` = "#D9D9D9")
marker_genes <- c("TRDC", "TRGC1", "TRGC2", "TRAC", "CD3E", "CD8A", "CD8B", "GZMK", "KLRF1", "FCGR3A", "NCAM1", "GNLY",
                  "KIT", "IL7R", "RORC", "ZBTB16", "KLRB1", "CX3CR1", "HSPA1A", "MKI67")
marker_genes <- marker_genes[marker_genes %in% rownames(sub)]

plot_data <- data.table(Embeddings(sub, "umap"))
plot_data[, is_tgd := sub$label %in% "T_gd" & sub$dataset != "s1"]   # published T_gd (non-S1)
plot_data <- plot_data[order(is_tgd)]                                 # draw T_gd cells on top

pdf(figure_file("umaps.pdf"), width = 20, height = 8)
# page 1: final annotation | subclusters at res 1.2
final_plot <- DimPlot(sub, group.by = "final_label", shuffle = TRUE, pt.size = 0.6, label = TRUE, repel = TRUE, cols = label_colors) +
              ggtitle(sprintf("Final annotation (%d cells; %s)", ncol(sub), Sys.getenv("TAG", "r2")))
res12_plot <- DimPlot(sub, group.by = "RNA_snn_res.1.2", pt.size = 0.6, label = TRUE) + ggtitle("subclusters, res 1.2")
print(final_plot | res12_plot)
# page 2: where are the published T_gd cells | dataset
tgd_plot <- ggplot(plot_data, aes(umap_1, umap_2, colour = is_tgd)) +
            geom_point(size = 0.6, stroke = 0) +
            coord_equal() +
            theme_classic() +
            scale_colour_manual(values = c(`FALSE` = "grey85", `TRUE` = "#A65628"), labels = c("other", "published T_gd"), name = NULL) +
            ggtitle("Where is T_gd?")
dataset_plot <- DimPlot(sub, group.by = "dataset", shuffle = TRUE, pt.size = 0.6) + ggtitle("dataset")
print(tgd_plot | dataset_plot)
# page 3: subclusters at res 0.4 | res 0.8
res04_plot <- DimPlot(sub, group.by = "RNA_snn_res.0.4", pt.size = 0.6, label = TRUE) + ggtitle("res 0.4")
res08_plot <- DimPlot(sub, group.by = "RNA_snn_res.0.8", pt.size = 0.6, label = TRUE) + ggtitle("res 0.8")
print(res04_plot | res08_plot)
dev.off()

pdf(figure_file("markers.pdf"), width = 18, height = 12)
print(FeaturePlot(sub, features = marker_genes, ncol = 5, pt.size = 0.1, order = TRUE) &
      theme(plot.title = element_text(size = 9), axis.text = element_blank()))
print(DotPlot(sub, features = marker_genes, group.by = "RNA_snn_res.1.2") + RotatedAxis() + ggtitle("markers by subcluster (res 1.2)"))
dev.off()
cat("done\n")
