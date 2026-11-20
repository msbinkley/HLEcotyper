# =============================================================================
# v2_C15_subumap_r14_cd8_nk_ilc_tgd.R
#
# Purpose (author request, 2026-10-25):
#   "make me a subumap of tgd nk ilc cd8 clusters, i want to see where tgd is"  (round 14 = the current standard)
#
# Cells:
#   - every cell of the round-14 clusters whose type is T_CD8, NK or ILC (Louvain res 0.6 clusters 4, 9 and 14;
#     12,763 cells). All 717 round-14 T_gd cells sit in clusters 4 (308) and 9 (409), so none is left out.
#
# Recipe: the same as the round-13 sub-UMAP v2_C11 (author's choice for T_gd sub-UMAPs, v2_C5c / v2_C5g):
#   LogNormalize -> 4,000 variable genes picked per dataset, ranked by how many datasets pick them
#   (author 2026-10-25: "i think i pick 4000 variable genes"; as the round-14 recluster v2_C13; v2_C11 used 3,000)
#   -> scale -> PCA 30 -> Harmony on dataset + chemistry + sample (theta 4 / 4 / 2)
#   -> Louvain at res 0.4, 0.8 and 1.2 -> UMAP.
#   (SCTransform is not used: S1 is one dataset = one chemistry = one sample and did not integrate with it.)
#
# Counts: the step-C count matrix (has TRDC / TRGC1 / TRGC2 / TRDV / TRGV genes, which the clustered object lacks).
#
# Markers (author rule): FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) on log-normalised RNA,
#   no downsampling; top 50 per subcluster ranked by p_adj, then avg_log2FC; report how many of the 50 have p_adj < 0.05.
#
# Nothing is removed or relabelled here (T_gd is protected). The author reads the plots.
#
# Run on Sherlock:  sbatch v2_run.sh v2_C15_subumap_r14_cd8_nk_ilc_tgd.R
# Output: output/subumap/2026-10-25_r14_sub_cd8_nk_ilc_tgd_*
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(harmony)
  library(ggplot2)
})
options(future.globals.maxSize = 16 * 1024^3)

work_dir <- "${PROJECT_ROOT}/cHL_v2"
out_dir  <- file.path(work_dir, "output/subumap")
out_file <- function(name) file.path(out_dir, paste0("2026-10-25_r14_sub_cd8_nk_ilc_tgd_", name))
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

seed        <- 42                       # same seed as the previous T_gd sub-UMAP
resolutions <- c("0.4", "0.8", "1.2")
harmony_theta <- c(dataset = 4, chemistry = 4, sample = 2)   # author choice B, 2026-10-16


# Step 1. Pick the cells: round-14 clusters 4 (T_CD8), 9 (NK), 14 (ILC) only
cells_r14 <- fread(cmd = sprintf("gzip -dc %s", file.path(work_dir, "output/cpA14_cells.tsv.gz")))
cells_r14[, whole_cluster := as.character(RNA_snn_res.0.6)]

cluster_type <- cells_r14[, .N, by = .(whole_cluster, final_label)]
cluster_type <- cluster_type[order(whole_cluster, -N)]
cluster_type <- cluster_type[, .(type = final_label[1]), by = whole_cluster]
wanted_clusters <- cluster_type[type %in% c("T_CD8", "NK", "ILC"), whole_cluster]
say("round-14 clusters whose type is T_CD8 / NK / ILC: %s", paste(wanted_clusters, collapse = ", "))
stopifnot(setequal(wanted_clusters, c("4", "9", "14")))

in_wanted_cluster <- cells_r14$whole_cluster %in% wanted_clusters
is_tgd            <- cells_r14$final_label == "T_gd"
chosen <- cells_r14[in_wanted_cluster]
say("cells chosen: %d; T_gd left out because they sit in other clusters: %d",
    nrow(chosen), sum(is_tgd & !in_wanted_cluster))
print(chosen[, .N, by = final_label][order(-N)])


# Step 2. Build the Seurat object from the step-C count matrix (gene symbols x cells)
counts_and_meta <- readRDS(file.path(work_dir, "output/stepC/cHL_v2_r10_counts_universe.rds"))
stopifnot(all(chosen$cell %in% colnames(counts_and_meta$counts)))
counts <- counts_and_meta$counts[, chosen$cell]
rm(counts_and_meta)
invisible(gc())

meta <- data.frame(
  dataset       = chosen$dataset,
  sample        = chosen$sample,
  chemistry     = chosen$chemistry,
  published     = chosen$label,
  final_label   = chosen$final_label,
  whole_cluster = chosen$whole_cluster,
  row.names     = chosen$cell
)
sub <- CreateSeuratObject(counts, meta.data = meta)
rm(counts)
invisible(gc())


# Step 3. Normalise
sub <- NormalizeData(sub, verbose = FALSE)


# Step 4. Variable genes: N_HVG per dataset, then rank genes by how many datasets pick them
#         (ties broken by the gene's median rank across datasets); keep the top N_HVG
N_HVG <- 4000                           # author 2026-10-25 (round-14 recluster also used 4,000)
hvg_per_dataset <- list()
for (ds in unique(sub$dataset)) {
  ds_cells  <- colnames(sub)[sub$dataset == ds]
  ds_object <- CreateSeuratObject(LayerData(sub, assay = "RNA", layer = "counts")[, ds_cells])
  ds_object <- FindVariableFeatures(ds_object, nfeatures = N_HVG, verbose = FALSE)
  hvg_per_dataset[[ds]] <- VariableFeatures(ds_object)
}
all_hvg <- unique(unlist(hvg_per_dataset))
n_datasets  <- numeric(length(all_hvg))
median_rank <- numeric(length(all_hvg))
for (i in seq_along(all_hvg)) {
  ranks <- c()
  for (ds in names(hvg_per_dataset)) {
    ranks <- c(ranks, match(all_hvg[i], hvg_per_dataset[[ds]]))
  }
  n_datasets[i]  <- sum(!is.na(ranks))
  median_rank[i] <- median(ranks, na.rm = TRUE)
}
hvg <- all_hvg[order(-n_datasets, median_rank)][1:N_HVG]
VariableFeatures(sub) <- hvg


# Step 5. Scale, PCA 30, Harmony (dataset + chemistry + sample), neighbours, clusters, UMAP
sub <- ScaleData(sub, features = hvg, verbose = FALSE)
sub <- RunPCA(sub, features = hvg, npcs = 30, seed.use = seed, verbose = FALSE)
set.seed(seed)
sub <- RunHarmony(sub, group.by.vars = names(harmony_theta), theta = unname(harmony_theta),
                  reduction.use = "pca", dims.use = 1:30, verbose = FALSE)
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (res in resolutions) {
  sub <- FindClusters(sub, resolution = as.numeric(res), random.seed = seed, verbose = FALSE)
}
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:30, seed.use = seed, verbose = FALSE)
saveRDS(sub, out_file("obj.rds"))
say("object saved: %d cells", ncol(sub))


# Step 6. Composition of every subcluster: cells, current labels, % T_gd, % of all T_gd, datasets,
#         and the % of cells expressing the gamma-delta TCR genes (TRDC, TRGC1, TRGC2)
tcr_genes <- intersect(c("TRDC", "TRGC1", "TRGC2", "TRAC"), rownames(sub))
expr <- LayerData(sub, assay = "RNA", layer = "data")
cell_table <- data.table(cell = colnames(sub), dataset = sub$dataset, final_label = sub$final_label)
for (g in tcr_genes) {
  cell_table[[paste0("expr_", g)]] <- as.numeric(expr[g, ] > 0)
}
total_tgd <- sum(cell_table$final_label == "T_gd")

composition <- list()
for (res in resolutions) {
  cell_table$subcluster <- as.character(sub[[paste0("RNA_snn_res.", res)]][, 1])
  for (sc in sort(unique(as.integer(cell_table$subcluster)))) {
    in_sc <- cell_table[subcluster == as.character(sc)]
    label_counts <- sort(table(in_sc$final_label), decreasing = TRUE)
    row <- data.table(
      res             = res,
      subcluster      = sc,
      cells           = nrow(in_sc),
      T_gd            = sum(in_sc$final_label == "T_gd"),
      pct_T_gd        = round(100 * mean(in_sc$final_label == "T_gd"), 1),
      pct_of_all_T_gd = round(100 * sum(in_sc$final_label == "T_gd") / total_tgd, 1),
      majority_label  = names(label_counts)[1],
      labels          = paste(names(label_counts), label_counts, sep = "=", collapse = " "),
      datasets        = paste(names(table(in_sc$dataset)), table(in_sc$dataset), sep = "=", collapse = " ")
    )
    for (g in tcr_genes) {
      row[[paste0("pct_", g, "_pos")]] <- round(100 * mean(in_sc[[paste0("expr_", g)]]), 1)
    }
    composition[[length(composition) + 1]] <- row
  }
}
composition <- rbindlist(composition)
fwrite(composition, out_file("composition.tsv"), sep = "\t")
print(composition[res == "0.8"])
print(composition[res == "1.2"])


# Step 7. Markers per subcluster (author rule), at res 0.8 and 1.2
DefaultAssay(sub) <- "RNA"
for (res in c("0.8", "1.2")) {
  Idents(sub) <- paste0("RNA_snn_res.", res)
  markers <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  top50 <- markers %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 50)
  fwrite(top50, out_file(sprintf("top50_markers_res%s.tsv", res)), sep = "\t")

  summary_rows <- list()
  for (sc in unique(top50$cluster)) {
    sc_markers <- top50[top50$cluster == sc, ]
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      subcluster = sc,
      n_top50    = nrow(sc_markers),
      n_padj_lt_0.05 = sum(sc_markers$p_val_adj < 0.05),
      top15      = paste(head(sc_markers$gene, 15), collapse = " ")
    )
  }
  fwrite(rbindlist(summary_rows), out_file(sprintf("top50_markers_res%s_summary.tsv", res)), sep = "\t")
  say("markers done at res %s", res)
}


# Step 8. Cell table (for plotting and annotation on the Mac): UMAP, subclusters, labels, TCR genes
umap <- Embeddings(sub, "umap")
cell_out <- data.table(cell = colnames(sub), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
                       dataset = sub$dataset, sample = sub$sample, published = sub$published,
                       final_label = sub$final_label, whole_cluster = sub$whole_cluster)
for (res in resolutions) {
  cell_out[[paste0("sub_res", res)]] <- as.character(sub[[paste0("RNA_snn_res.", res)]][, 1])
}
fwrite(cell_out, out_file("cells.tsv.gz"), sep = "\t")


# Step 9. Figures
marker_genes <- c("TRDC", "TRGC1", "TRGC2", "TRDV1", "TRDV2", "TRGV9", "TRAC", "CD3E", "CD8A", "CD8B", "GZMK", "GZMB",
                  "KLRF1", "FCGR3A", "NCAM1", "GNLY", "NKG7", "KIT", "IL7R", "RORC", "ZBTB16", "KLRB1", "CX3CR1", "MKI67")
marker_genes <- marker_genes[marker_genes %in% rownames(sub)]
label_colors <- c(ILC = "#00C19A", NK = "#FDBF6F", T_CD8 = "#999999", T_gd = "#A65628", T_CD4 = "#B2DF8A")

plot_data <- data.table(Embeddings(sub, "umap"), keep.rownames = "cell")
plot_data$is_tgd <- sub$final_label == "T_gd"
plot_data <- plot_data[order(is_tgd)]                    # draw T_gd cells on top

pdf(out_file("umaps.pdf"), width = 12, height = 9)
print(DimPlot(sub, group.by = "final_label", shuffle = TRUE, pt.size = 0.3, label = TRUE, repel = TRUE, cols = label_colors) +
      ggtitle(sprintf("Sub-UMAP of round-14 clusters 4 / 9 / 14: %d cells; round-14 final labels", ncol(sub))))
print(ggplot(plot_data, aes(umap_1, umap_2, colour = is_tgd)) +
      geom_point(size = 0.3, stroke = 0) +
      scale_colour_manual(values = c(`FALSE` = "grey85", `TRUE` = "#A65628"), labels = c("other", "current T_gd"), name = NULL) +
      coord_equal() + theme_classic() +
      ggtitle(sprintf("Where is T_gd? (%d current T_gd cells)", sum(plot_data$is_tgd))))
for (res in resolutions) {
  print(DimPlot(sub, group.by = paste0("RNA_snn_res.", res), label = TRUE, pt.size = 0.3) + ggtitle(paste("subclusters, res", res)))
}
print(DimPlot(sub, group.by = "dataset", shuffle = TRUE, pt.size = 0.3) + ggtitle("dataset"))
print(FeaturePlot(sub, features = marker_genes, ncol = 6, pt.size = 0.1, order = TRUE) &
      theme(plot.title = element_text(size = 9), axis.text = element_blank()))
for (res in c("0.8", "1.2")) {
  print(DotPlot(sub, features = marker_genes, group.by = paste0("RNA_snn_res.", res)) + RotatedAxis() +
        ggtitle(paste("markers by subcluster, res", res)))
}
print(DotPlot(sub, features = marker_genes, group.by = "final_label") + RotatedAxis() + ggtitle("markers by current final label"))
dev.off()
say("done")
