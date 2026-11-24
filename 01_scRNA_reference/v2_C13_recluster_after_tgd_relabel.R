# =============================================================================
# v2_C13_recluster_after_tgd_relabel.R
#
# Purpose (author, 2026-10-23):
#   "If you removed 11 T_gd cells from the T_CD4 cluster, you should recluster the UMAP."
#   Take the round-13 object, drop the 11 removed T_gd cells, attach the new labels from the T_gd sub-UMAP
#   re-annotation (v2_C12_relabel_from_tgd_subumap.R), and recluster + re-run the UMAP with the SAME recipe as the
#   round-13 final recluster (v2_C8_round13_loop.R, step 3):
#     4,000 variable genes picked per dataset and ranked by how many datasets pick them -> scale -> PCA 50
#     (author 2026-10-23: 4,000 variable genes, not round 13's 3,000)
#     -> Harmony on dataset + chemistry + sample (theta 4 / 4 / 2) -> neighbours on 30 Harmony dims
#     -> Louvain at res 0.3, 0.6, 1.0 (seed 1399, as round 13) -> UMAP (seed 42), integrated and unintegrated
#     -> iLISI (inverse Simpson index of dataset / chemistry over 30 neighbours), as round 13.
#   Labels are NOT changed here; the script only reports how the new clusters line up with the new labels.
#
# Inputs (OAK):  output/cHL_merged_harmony_v2_r13.rds
#                output/stepC/inputs/2026-10-23_round13_labels_tgd_subumap.tsv.gz  (112,364 cells)
# Outputs (OAK): output/cHL_merged_harmony_v2_r13b.rds
#                output/cpA13b_cells.tsv.gz  (same columns as cpA13_cells.tsv.gz: UMAP, clusters, final_label, QC)
#                output/cpA13b_cluster_by_final_label.tsv, output/cpA13b_hvg_4000.tsv
# Run on Sherlock:  sbatch v2_run.sh v2_C13_recluster_after_tgd_relabel.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
})
work_dir <- "${PROJECT_ROOT}/cHL_v2"
out_dir  <- file.path(work_dir, "output")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

cluster_seed  <- 1399                                        # same seed as the round-13 final recluster
umap_seed     <- 42
n_hvg         <- 4000                                        # author 2026-10-23 (round 13 used 3,000)
harmony_theta <- c(dataset = 4, chemistry = 4, sample = 2)   # same as round 13
gene_table <- fread(file.path(work_dir, "ref/gene_table_v47.tsv"))


# Step 1. Read the round-13 object and the new labels; keep only the cells that have a new label
obj <- readRDS(file.path(out_dir, "cHL_merged_harmony_v2_r13.rds"))
new_labels <- fread(cmd = sprintf("gzip -dc %s", file.path(out_dir, "stepC/inputs/2026-10-23_round13_labels_tgd_subumap.tsv.gz")))
removed_cells <- setdiff(colnames(obj), new_labels$cell)
say("round 13: %d cells; new label file: %d cells; removed: %d", ncol(obj), nrow(new_labels), length(removed_cells))
stopifnot(length(removed_cells) == 11)
stopifnot(all(obj$r13_label[removed_cells] == "T_gd"))

obj <- subset(obj, cells = new_labels$cell)
obj$r13_cluster_res0.6 <- as.character(obj$RNA_snn_res.0.6)     # keep the round-13 clusters for comparison
obj$final_label <- new_labels$final_label[match(colnames(obj), new_labels$cell)]


# Step 2. Variable genes: 4,000 per dataset, then rank genes by how many datasets pick them
#         (ties broken by the median rank across datasets); keep the top 4,000
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
all_hvg <- names(table(unlist(hvg_per_dataset)))                  # same gene order as round 13 (table() sorts names)
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
fwrite(data.table(ensg = hvg, symbol = gene_table$symbol[match(hvg, gene_table$ensg)]), file.path(out_dir, "cpA13b_hvg_4000.tsv"), sep = "\t")


# Step 3. Scale, PCA 50, Harmony, neighbours, clusters at res 0.3 / 0.6 / 1.0
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = cluster_seed, verbose = FALSE)
set.seed(cluster_seed)
obj <- RunHarmony(obj, group.by.vars = names(harmony_theta), theta = unname(harmony_theta),
                  reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (res in c(0.3, 0.6, 1.0)) {
  obj <- FindClusters(obj, resolution = res, random.seed = cluster_seed, verbose = FALSE)
}


# Step 4. UMAPs (unintegrated on PCA, integrated on Harmony), seed 42 as round 13
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = umap_seed, reduction.name = "umap_unintegrated", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = umap_seed, verbose = FALSE)


# Step 5. iLISI: for each cell, inverse Simpson index of dataset (and chemistry) among its 30 nearest neighbours
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
say("round 13b: %d cells, %d clusters at res 0.6; median iLISI dataset %.2f", ncol(obj),
    length(unique(obj$RNA_snn_res.0.6)), median(obj$ilisi_dataset))


# Step 6. Per-cell table (same columns as cpA13_cells.tsv.gz, final_label = new label)
umap <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
keep_columns <- c("dataset", "sample", "chemistry", "label", "s1_cell_type", "ref_exclude", "r12_cluster", "r12_label", "final_label",
                  "RNA_snn_res.0.3", "RNA_snn_res.0.6", "RNA_snn_res.1", "ilisi_dataset", "ilisi_chemistry",
                  "nCount_RNA", "nFeature_RNA", "nCount_orig", "nFeature_orig", "percent.mt", "r13_cluster_res0.6")
keep_columns <- intersect(keep_columns, names(obj@meta.data))
cells <- data.table(obj@meta.data[, keep_columns], keep.rownames = "cell")
cells$UMAP_1  <- umap[, 1]
cells$UMAP_2  <- umap[, 2]
cells$UMAPu_1 <- umap_unintegrated[, 1]
cells$UMAPu_2 <- umap_unintegrated[, 2]
fwrite(cells, file.path(out_dir, "cpA13b_cells.tsv.gz"), sep = "\t")


# Step 7. How the new clusters (res 0.6) line up with the new labels
by_cluster <- cells[, .N, by = .(cluster = RNA_snn_res.0.6, final_label)]
by_cluster[, pct := round(100 * N / sum(N), 1), by = cluster]
by_cluster <- by_cluster[order(cluster, -N)]
fwrite(by_cluster, file.path(out_dir, "cpA13b_cluster_by_final_label.tsv"), sep = "\t")
print(by_cluster[, .SD[1], by = cluster][order(pct)])

old_vs_new <- cells[, .N, by = .(r13_cluster_res0.6, new_cluster = RNA_snn_res.0.6)]
print(old_vs_new[order(r13_cluster_res0.6, -N)][, .SD[1], by = r13_cluster_res0.6])

saveRDS(obj, file.path(out_dir, "cHL_merged_harmony_v2_r13b.rds"))
say("saved round 13b")
