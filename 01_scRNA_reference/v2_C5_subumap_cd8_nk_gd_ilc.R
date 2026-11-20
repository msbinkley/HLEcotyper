# =============================================================================
# v2_C5_subumap_cd8_nk_gd_ilc.R
#
# Purpose (author request, 2026-10-14):
#   "make me a sub-UMAP on T_CD8, T_gd, NK and ILC; I need to know where T_gd is"
#   (T_gd forms no cluster on the whole-object UMAP at res 0.6 - 4).
#   Sub-UMAP of T_CD8 / T_gd / NK / ILC on the final rolling object (r11_rolling).
#
# Cells:
#   - every cell of r11_rolling clusters 4 (T_CD8), 8 (NK) and 16 (ILC)
#   - plus every cell published as T_CD8 / T_gd / NK / ILC, wherever it sits
#
# Counts: the step-C count matrix (17,823 genes incl. TRDC / TRGC1 / TRGC2 / TRAC, which the clustered object lacks).
#
# Recipe = v2_C1_subumap.R:
#   SCTransform v2 (glmGamPoi) -> PCA 30 -> Harmony (dataset + chemistry + sample, seed 42)
#   -> Louvain res 0.4 / 0.8 / 1.2 -> UMAP.
#   Markers on log-normalised RNA (Wilcoxon), res 0.8.
#   (This is v1. S1 stayed a separate island with SCTransform; v2_C5b redoes it with LogNormalize.)
#
# Run on Sherlock.
# Output: output/subumap/r11rolling_sub_CD8_NK_gd_ILC_*  (object, composition, markers, PDF)
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})
options(future.globals.maxSize = 16 * 1024^3)

work_dir <- "${PROJECT_ROOT}/cHL_v2"
out_dir  <- file.path(work_dir, "output/subumap")
tag      <- "r11rolling_sub_CD8_NK_gd_ILC"
out_file <- function(name) file.path(out_dir, paste0(tag, "_", name))           # output/subumap/<tag>_<name>
say      <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")  # time-stamped log line


# Step 1. Pick the cells: whole-object clusters 4 / 8 / 16 + every cell published as T_CD8 / T_gd / NK / ILC
whole_object <- readRDS(file.path(work_dir, "output/cHL_merged_harmony_v2_r11_rolling.rds"))
cell_table <- data.table(whole_object@meta.data, keep.rownames = "cell")
cell_table <- cell_table[, .(cell, dataset, sample, chemistry, label, s1_cell_type)]
cell_table[, whole_cluster := as.character(whole_object$RNA_snn_res.0.6)]     # whole-object Louvain res 0.6 cluster
rm(whole_object)
invisible(gc())

in_wanted_cluster  <- cell_table$whole_cluster %in% c("4", "8", "16")
has_wanted_label   <- cell_table$label %in% c("T_CD8", "T_gd", "NK", "ILC")
chosen_cells <- cell_table$cell[in_wanted_cluster | has_wanted_label]


# Step 2. Build the Seurat object from the step-C count matrix
counts_and_meta <- readRDS(file.path(work_dir, "output/stepC/cHL_v2_r10_counts_universe.rds"))   # list(counts = genes(symbols) x cells, meta)
stopifnot(all(chosen_cells %in% colnames(counts_and_meta$counts)))

chosen_meta <- cell_table[match(chosen_cells, cell)]
# Label shown on plots: S1 has no published labels here; other cells with no label are "unlabelled"
chosen_meta[, label_show := fifelse(dataset == "s1", "S1 (unlabelled here)",
                                    fifelse(is.na(label), "unlabelled", label))]

sub <- CreateSeuratObject(counts_and_meta$counts[, chosen_cells],
                          meta.data = data.frame(chosen_meta, row.names = chosen_meta$cell))
rm(counts_and_meta)
invisible(gc())
label_counts <- table(sub$label_show)
say("%d cells: %s", ncol(sub), paste(names(label_counts), label_counts, sep = "=", collapse = " "))


# Step 3. SCTransform v2 -> PCA 30 -> Harmony -> neighbours -> Louvain at 3 resolutions -> UMAP (all seeds 42)
sub <- SCTransform(sub, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
set.seed(42)
sub <- RunHarmony(sub, group.by.vars = c("dataset", "chemistry", "sample"),   # default theta (2)
                  reduction.use = "pca", reduction.save = "harmony", verbose = FALSE)
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (res in c(0.4, 0.8, 1.2)) {
  sub <- FindClusters(sub, resolution = res, random.seed = 42, verbose = FALSE)   # adds SCT_snn_res.<res>
}
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)
saveRDS(sub, out_file("obj.rds"))


# Step 4. Composition: where are the T_gd cells?
#   One row per subcluster and resolution: cells, published T_gd count and %, % of all T_gd at that resolution,
#   majority published label (S1 cells ignored), top 3 labels with their %, cells per dataset.

# Helper: the composition row for the cells of one subcluster
summarise_subcluster <- function(label_show, dataset) {
  sorted_labels <- sort(table(label_show), decreasing = TRUE)
  n_top         <- min(3, uniqueN(label_show))                   # top 3 labels (fewer if the subcluster has fewer)
  labelled      <- label_show[label_show != "S1 (unlabelled here)"]
  list(
    cells    = length(label_show),
    T_gd     = sum(label_show == "T_gd"),
    pct_T_gd = round(100 * mean(label_show == "T_gd"), 1),
    majority = names(which.max(table(labelled))),
    labels   = paste(names(sorted_labels)[1:n_top],
                     round(100 * sorted_labels[1:n_top] / length(label_show)), collapse = "; "),
    datasets = paste(names(table(dataset)), table(dataset), sep = "=", collapse = " ")
  )
}

cell_meta <- data.table(sub@meta.data, keep.rownames = "cell")
composition <- list()
for (res_value in c("0.4", "0.8", "1.2")) {
  cell_meta[, subcluster := get(paste0("SCT_snn_res.", res_value))]
  res_table <- cell_meta[, summarise_subcluster(label_show, dataset), by = subcluster]   # subclusters in order of first appearance
  res_table[, res := res_value]
  res_table[, pct_of_all_T_gd := round(100 * T_gd / sum(T_gd), 1)]
  res_table <- res_table[order(-pct_T_gd)]                            # most T_gd-rich subcluster first
  composition[[length(composition) + 1]] <- res_table
}
composition <- rbindlist(composition)
setcolorder(composition, c("res", "subcluster", "cells", "T_gd", "pct_T_gd", "pct_of_all_T_gd", "majority"))
fwrite(composition, out_file("composition.tsv"), sep = "\t")
print(composition[res == "0.8"])
print(composition[res == "1.2"][1:6])


# Step 5. Markers per subcluster at res 0.8, on log-normalised RNA
#   (Wilcoxon; min.pct 0.1, at most 500 cells per subcluster (seed 42); keep p_adj < 0.05, top 30 by log2FC)
DefaultAssay(sub) <- "RNA"
sub <- NormalizeData(sub, verbose = FALSE)
Idents(sub) <- "SCT_snn_res.0.8"
markers <- data.table(FindAllMarkers(sub, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1,
                                     max.cells.per.ident = 500, random.seed = 42, verbose = FALSE))
significant <- markers[p_val_adj < 0.05]
significant <- significant[order(cluster, -avg_log2FC)]
top30 <- significant[, head(.SD, 30), by = cluster]                  # first 30 rows of each subcluster
fwrite(top30, out_file("top30_markers_res0.8.tsv"), sep = "\t")


# Step 6. Figures
marker_genes <- c("TRDC", "TRGC1", "TRGC2", "TRAC", "CD3E", "CD8A", "CD8B", "GZMK", "KLRF1", "FCGR3A", "NCAM1", "GNLY",
                  "KIT", "IL7R", "RORC", "ZBTB16", "KLRB1", "CX3CR1")
marker_genes <- marker_genes[marker_genes %in% rownames(sub)]
label_colors <- c(ILC = "#00C19A", NK = "#FDBF6F", T_CD8 = "#999999", T_gd = "#A65628", T_CD4 = "#B2DF8A",
                  T_Treg = "#A6CEE3", T_MKI67 = "#FADD4B", `S1 (unlabelled here)` = "#D9D9D9", unlabelled = "#D9D9D9")

plot_data <- data.table(Embeddings(sub, "umap"), keep.rownames = "cell")
plot_data[, is_tgd := sub$label_show == "T_gd"]
plot_data <- plot_data[order(is_tgd)]                                 # draw T_gd cells on top

pdf(out_file("umaps.pdf"), width = 12, height = 9)
print(DimPlot(sub, group.by = "label_show", shuffle = TRUE, pt.size = 0.3, label = TRUE, repel = TRUE, cols = label_colors) +
      ggtitle(sprintf("T_CD8 / T_gd / NK / ILC sub-UMAP: %d cells; published labels (S1 grey)", ncol(sub))))
print(ggplot(plot_data, aes(umap_1, umap_2, colour = is_tgd)) +
      geom_point(size = 0.3, stroke = 0) +
      scale_colour_manual(values = c(`FALSE` = "grey85", `TRUE` = "#A65628"), labels = c("other", "published T_gd"), name = NULL) +
      coord_equal() +
      theme_classic() +
      ggtitle(sprintf("Where is T_gd? (%d published T_gd cells)", sum(plot_data$is_tgd))))
for (res in c("0.4", "0.8", "1.2")) {
  print(DimPlot(sub, group.by = paste0("SCT_snn_res.", res), label = TRUE, pt.size = 0.3) + ggtitle(paste("subclusters, res", res)))
}
print(DimPlot(sub, group.by = "dataset", shuffle = TRUE, pt.size = 0.3) + ggtitle("dataset"))
print(FeaturePlot(sub, features = marker_genes, ncol = 6, pt.size = 0.1, order = TRUE) &
      theme(plot.title = element_text(size = 9), axis.text = element_blank()))
print(DotPlot(sub, features = marker_genes, group.by = "SCT_snn_res.0.8") + RotatedAxis() + ggtitle("markers by subcluster (res 0.8)"))
print(DotPlot(sub, features = marker_genes, group.by = "label_show") + RotatedAxis() + ggtitle("markers by published label"))
dev.off()
say("done")
