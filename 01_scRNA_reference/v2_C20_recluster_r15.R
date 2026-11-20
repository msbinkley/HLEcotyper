# =============================================================================
# v2_C20_recluster_r15.R
#
# Purpose (author, 2026-10-25):
#   "record c7 in subclusters, res 0.4 as tgd. record cells that were removed in subumap, then remove these cells from
#    r14 to r15. record cells in c7 to tgd in r15. recluster r15"
#   Round 15 = round 14 minus the 569 cells removed in the CD8 / NK / ILC / T_gd sub-UMAP chain (r14 -> r14b -> r14c);
#   T_gd = the 357 cells of res-0.4 c7 of the r14c sub-UMAP; all other cells keep their round-14 label.
#   The records and the label file are written on the Mac by v2_C19_record_subumap_removals_r15.R.
#
# Recluster: the SAME recipe and seed as the round-14 final recluster (v2_C14_unstable_loop_r13b.R, step 3; as the
#   round-13b recluster kept round 13's seed), so the only change is the cells and labels:
#     4,000 variable genes per dataset, ranked by how many datasets pick them -> scale -> PCA 50
#     -> Harmony on dataset + chemistry + sample (theta 4 / 4 / 2) -> neighbours on 30 Harmony dims
#     -> Louvain at res 0.3 / 0.6 / 1.0 (seed 1499, as round 14) -> UMAP (seed 42), integrated and unintegrated
#     -> iLISI (inverse Simpson index of dataset / chemistry over 30 neighbours).
#   Labels are NOT changed by the recluster; the script only reports how the new clusters line up with them.
#   No unstable-cell loop is run here (not asked).
#
# Inputs (OAK):  output/cHL_merged_harmony_v2_r14.rds
#                output/stepC/inputs/2026-10-25_round15_labels_pre.tsv.gz (v2_C19: cell, dataset, r14_label, final_label)
# Outputs (OAK): output/cHL_merged_harmony_v2_r15.rds
#                output/2026-10-25_cpA15_cells.tsv.gz            (same columns as cpA14_cells.tsv.gz + r14 cluster / label)
#                output/2026-10-25_cpA15_cluster_by_final_label.tsv, output/2026-10-25_cpA15_hvg_4000.tsv
#                output/2026-10-25_cpA15_umaps.pdf
#                output/stepC/inputs/2026-10-25_round15_labels.tsv.gz (same layout as the round-14 label file)
# Run on Sherlock:  sbatch v2_run.sh v2_C20_recluster_r15.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
  library(ggplot2)
})
work_dir <- "${PROJECT_ROOT}/cHL_v2"
out_dir  <- file.path(work_dir, "output")
out_file <- function(name) file.path(out_dir, paste0("2026-10-25_cpA15_", name))
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

cluster_seed   <- 1499                                       # same seed as the round-14 final recluster
umap_seed      <- 42                                         # same as round 14
n_hvg          <- 4000                                       # same as round 14
harmony_theta  <- c(dataset = 4, chemistry = 4, sample = 2)  # same as round 14
gene_table <- fread(file.path(work_dir, "ref/gene_table_v47.tsv"))


# Step 1. Read the round-14 object and the round-15 labels (v2_C19)
obj <- readRDS(file.path(out_dir, "cHL_merged_harmony_v2_r14.rds"))
r15_labels <- fread(cmd = sprintf("gzip -dc %s", file.path(out_dir, "stepC/inputs/2026-10-25_round15_labels_pre.tsv.gz")))
stopifnot(all(r15_labels$cell %in% colnames(obj)))
stopifnot(all(obj$final_label[r15_labels$cell] == r15_labels$r14_label))     # same round-14 labels as on Box
say("round 14: %d cells; round 15: %d cells; removed: %d", ncol(obj), nrow(r15_labels), ncol(obj) - nrow(r15_labels))


# Step 2. Old and new label of every round-14 cell (removed cells = "REMOVED")
old_label <- obj$final_label
names(old_label) <- colnames(obj)
new_label <- rep("REMOVED", length(old_label))
names(new_label) <- names(old_label)
new_label[r15_labels$cell] <- r15_labels$final_label
changes <- data.table(old_label = old_label, new_label = new_label)[old_label != new_label, .N, by = .(old_label, new_label)]
print(changes[order(-N)])


# Step 3. Drop the removed cells and attach the new labels; keep the round-14 clusters for comparison
kept_cells <- names(new_label)[new_label != "REMOVED"]
obj <- subset(obj, cells = kept_cells)
obj$r14_cluster_res0.6 <- as.character(obj$RNA_snn_res.0.6)
obj$r14_label <- old_label[colnames(obj)]
obj$final_label <- new_label[colnames(obj)]
say("round 15: %d cells; T_gd %d", ncol(obj), sum(obj$final_label == "T_gd"))


# Step 4. Variable genes: 4,000 per dataset, then rank genes by how many datasets pick them
#         (ties broken by the median rank across datasets); keep the top 4,000 (as v2_C14 recluster())
obj@reductions <- list()
obj@graphs <- list()
obj@neighbors <- list()
hvg_per_dataset <- list()
for (ds in unique(obj$dataset)) {
  ds_cells  <- colnames(obj)[obj$dataset == ds]
  ds_object <- CreateSeuratObject(LayerData(obj, assay = "RNA", layer = "counts")[, ds_cells])
  ds_object <- FindVariableFeatures(ds_object, nfeatures = n_hvg, verbose = FALSE)
  hvg_per_dataset[[ds]] <- VariableFeatures(ds_object)
}
all_hvg <- names(table(unlist(hvg_per_dataset)))                  # same gene order as round 14 (table() sorts names)
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
hvg <- all_hvg[order(-n_datasets, median_rank)][1:n_hvg]
DefaultAssay(obj) <- "RNA"
VariableFeatures(obj) <- hvg
fwrite(data.table(ensg = hvg, symbol = gene_table$symbol[match(hvg, gene_table$ensg)]), out_file("hvg_4000.tsv"), sep = "\t")


# Step 5. Scale, PCA 50, Harmony, neighbours, clusters at res 0.3 / 0.6 / 1.0
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = cluster_seed, verbose = FALSE)
set.seed(cluster_seed)
obj <- RunHarmony(obj, group.by.vars = names(harmony_theta), theta = unname(harmony_theta),
                  reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (res in c(0.3, 0.6, 1.0)) {
  obj <- FindClusters(obj, resolution = res, random.seed = cluster_seed, verbose = FALSE)
}


# Step 6. UMAPs (unintegrated on PCA, integrated on Harmony), seed 42 as round 14
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = umap_seed, reduction.name = "umap_unintegrated", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = umap_seed, verbose = FALSE)


# Step 7. iLISI: for each cell, inverse Simpson index of dataset (and chemistry) among its 30 nearest neighbours
neighbours <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, k.param = 30, return.neighbor = TRUE, verbose = FALSE)@neighbors[[1]]
neighbour_index <- Indices(neighbours)
inverse_simpson <- function(labels) {
  apply(neighbour_index, 1, function(row) {
    shares <- table(labels[row]) / length(row)
    1 / sum(shares^2)
  })
}
obj$ilisi_dataset   <- inverse_simpson(obj$dataset)
obj$ilisi_chemistry <- inverse_simpson(obj$chemistry)
say("round 15: %d cells, %d clusters at res 0.6; median iLISI dataset %.2f; T_gd %d", ncol(obj),
    length(unique(obj$RNA_snn_res.0.6)), median(obj$ilisi_dataset), sum(obj$final_label == "T_gd"))


# Step 8. Per-cell table (same columns as cpA14_cells.tsv.gz, plus the round-14 cluster and label)
umap <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
keep_columns <- c("dataset", "sample", "chemistry", "label", "s1_cell_type", "ref_exclude", "r12_cluster", "r12_label", "final_label",
                  "RNA_snn_res.0.3", "RNA_snn_res.0.6", "RNA_snn_res.1", "ilisi_dataset", "ilisi_chemistry",
                  "nCount_RNA", "nFeature_RNA", "nCount_orig", "nFeature_orig", "percent.mt",
                  "r14_cluster_res0.6", "r14_label")
keep_columns <- intersect(keep_columns, names(obj@meta.data))
cells <- data.table(obj@meta.data[, keep_columns], keep.rownames = "cell")
cells$UMAP_1  <- umap[, 1]
cells$UMAP_2  <- umap[, 2]
cells$UMAPu_1 <- umap_unintegrated[, 1]
cells$UMAPu_2 <- umap_unintegrated[, 2]
fwrite(cells, out_file("cells.tsv.gz"), sep = "\t")


# Step 9. How the new clusters (res 0.6) line up with the labels, and with the round-14 clusters
by_cluster <- cells[, .N, by = .(cluster = RNA_snn_res.0.6, final_label)]
by_cluster[, pct := round(100 * N / sum(N), 1), by = cluster]
by_cluster <- by_cluster[order(cluster, -N)]
fwrite(by_cluster, out_file("cluster_by_final_label.tsv"), sep = "\t")
print(by_cluster[, .SD[1], by = cluster][order(pct)])
print(by_cluster[final_label == "T_gd"][order(-N)])

old_vs_new <- cells[, .N, by = .(r14_cluster_res0.6, new_cluster = RNA_snn_res.0.6)]
print(old_vs_new[order(r14_cluster_res0.6, -N)][, .SD[1], by = r14_cluster_res0.6])


# Step 10. Label file for step C, in the same layout as 2026-10-23_round14_labels.tsv.gz
majority <- by_cluster[, .(cluster_type = final_label[1]), by = cluster]
label_file <- cells[, .(cell, dataset, published_label = label, cluster = RNA_snn_res.0.6, final_label)]
label_file$cluster_type <- majority$cluster_type[match(label_file$cluster, majority$cluster)]
setcolorder(label_file, c("cell", "dataset", "published_label", "cluster", "cluster_type", "final_label"))
label_path <- file.path(out_dir, "stepC/inputs/2026-10-25_round15_labels.tsv")
fwrite(label_file, label_path, sep = "\t")
system(sprintf("gzip -f %s", label_path))


# Step 11. Review UMAPs: labels, clusters at res 0.6 / 1.0, where T_gd is, dataset
plot_data <- data.table(cell = colnames(obj), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2], is_tgd = obj$final_label == "T_gd")
plot_data <- plot_data[order(is_tgd)]                    # draw T_gd cells on top
pdf(out_file("umaps.pdf"), width = 12, height = 9)
print(DimPlot(obj, group.by = "final_label", shuffle = TRUE, raster = TRUE, label = TRUE, repel = TRUE) +
      ggtitle(sprintf("Round 15: %d cells; final labels (T_gd = r14c sub-UMAP res-0.4 c7)", ncol(obj))))
print(DimPlot(obj, group.by = "RNA_snn_res.0.6", raster = TRUE, label = TRUE) + ggtitle("Round 15 clusters, res 0.6"))
print(DimPlot(obj, group.by = "RNA_snn_res.1", raster = TRUE, label = TRUE) + ggtitle("Round 15 clusters, res 1.0"))
print(ggplot(plot_data, aes(UMAP_1, UMAP_2, colour = is_tgd)) +
      geom_point(size = 0.2, stroke = 0) +
      scale_colour_manual(values = c(`FALSE` = "grey85", `TRUE` = "#A65628"), labels = c("other", "T_gd"), name = NULL) +
      coord_equal() + theme_classic() +
      ggtitle(sprintf("Where is T_gd? (%d cells)", sum(plot_data$is_tgd))))
print(DimPlot(obj, group.by = "dataset", shuffle = TRUE, raster = TRUE) + ggtitle("dataset"))
dev.off()

saveRDS(obj, file.path(out_dir, "cHL_merged_harmony_v2_r15.rds"))
say("saved round 15")
