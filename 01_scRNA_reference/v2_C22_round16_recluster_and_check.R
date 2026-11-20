# =============================================================================
# v2_C22_round16_recluster_and_check.R
#
# Purpose (author, 2026-10-26):
#   "Remove the 162 cells and the 7 cells sit outside of HRS cluster. Then recluster and reidentify unstable cells"
#   The 162 unstable cells and the 7 HRS cells come from the round-15 check (v2_C21, seed 1501,
#   output/2026-10-26_cpA15check_unstable_cells.tsv.gz: unstable == TRUE, and label == HRS outside an HRS cluster).
#   The author removes the 7 HRS cells although HRS is protected (author decision, 2026-10-26).
#
# Part A. Round 16 = round 15 minus these 169 cells; labels unchanged. Recluster with the SAME recipe and seed as the
#   round-14 / 15 final recluster:
#     4,000 variable genes per dataset, ranked by how many datasets pick them -> scale -> PCA 50
#     -> Harmony on dataset + chemistry + sample (theta 4 / 4 / 2) -> neighbours on 30 Harmony dims
#     -> Louvain at res 0.3 / 0.6 / 1.0 (seed 1499) -> UMAP (seed 42), integrated and unintegrated -> iLISI.
# Part B. One unstable-cell round on round 16 (report only, as v2_C21) with a NEW seed (1502):
#   recluster at res 0.6; cluster type = most common round-16 label; unstable = label != cluster type;
#   T_gd and HRS protected (mismatches reported separately). Nothing is removed; the author decides.
#
# Inputs (OAK):  output/cHL_merged_harmony_v2_r15.rds, output/2026-10-26_cpA15check_unstable_cells.tsv.gz
# Outputs (OAK): output/2026-10-26_cpA16_removed_cells.tsv.gz       (the 169 cells, with reason)
#                output/cHL_merged_harmony_v2_r16.rds
#                output/2026-10-26_cpA16_cells.tsv.gz, _cluster_by_final_label.tsv, _hvg_4000.tsv
#                output/stepC/inputs/2026-10-26_round16_labels.tsv.gz
#                output/2026-10-26_cpA16check_{summary,by_type,moves,by_dataset,cluster_types}.tsv, _unstable_cells.tsv.gz
# Run on Sherlock:  sbatch v2_run.sh v2_C22_round16_recluster_and_check.R
# UMAP figures are drawn on the Mac from the cell table (DimPlot on this object defaults to the unintegrated UMAP).
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
out_file   <- function(name) file.path(out_dir, paste0("2026-10-26_cpA16_", name))
check_file <- function(name) file.path(out_dir, paste0("2026-10-26_cpA16check_", name))
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

cluster_seed   <- 1499                                       # same seed as the round-14 / 15 final recluster
check_seed     <- 1502                                       # new seed for the unstable check (1501 = round-15 check)
umap_seed      <- 42                                         # same as rounds 14 / 15
n_hvg          <- 4000                                       # same as rounds 14 / 15
harmony_theta  <- c(dataset = 4, chemistry = 4, sample = 2)  # same as rounds 14 / 15
protected      <- c("T_gd", "HRS")                           # rare types: not counted as unstable in part B
gene_table <- fread(file.path(work_dir, "ref/gene_table_v47.tsv"))


# ---------------------------------------------------------------------------------------------------------------------
# Part A. Round 16
# ---------------------------------------------------------------------------------------------------------------------

# Step 1. Read round 15 and the round-15 check; pick the 162 unstable cells and the 7 HRS cells outside an HRS cluster
obj <- readRDS(file.path(out_dir, "cHL_merged_harmony_v2_r15.rds"))
flagged <- fread(cmd = sprintf("gzip -dc %s", file.path(out_dir, "2026-10-26_cpA15check_unstable_cells.tsv.gz")))
unstable_cells <- flagged[unstable == TRUE]
hrs_outside    <- flagged[label == "HRS" & cluster_type != "HRS"]
stopifnot(nrow(unstable_cells) == 162)                     # numbers reported to the author
stopifnot(nrow(hrs_outside) == 7)
unstable_cells$reason <- "unstable in the round-15 check (seed 1501)"
hrs_outside$reason    <- "HRS outside the HRS cluster in the round-15 check (seed 1501)"
removed <- rbind(unstable_cells, hrs_outside)
stopifnot(!anyDuplicated(removed$cell), all(removed$cell %in% colnames(obj)))
fwrite(removed, out_file("removed_cells.tsv.gz"), sep = "\t")
print(removed[, .N, by = .(reason, label, cluster_type)][order(-N)])


# Step 2. Drop them; keep the round-15 clusters and labels for comparison (labels are not changed)
kept_cells <- setdiff(colnames(obj), removed$cell)
obj <- subset(obj, cells = kept_cells)
obj$r15_cluster_res0.6 <- as.character(obj$RNA_snn_res.0.6)
obj$r15_label <- obj$final_label
say("round 16: %d cells (removed %d); T_gd %d, HRS %d", ncol(obj), nrow(removed),
    sum(obj$final_label == "T_gd"), sum(obj$final_label == "HRS"))


# Step 3. Variable genes: 4,000 per dataset, then rank genes by how many datasets pick them
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
all_hvg <- names(table(unlist(hvg_per_dataset)))                  # same gene order as rounds 14 / 15 (table() sorts names)
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


# Step 4. Scale, PCA 50, Harmony, neighbours, clusters at res 0.3 / 0.6 / 1.0
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = cluster_seed, verbose = FALSE)
set.seed(cluster_seed)
obj <- RunHarmony(obj, group.by.vars = names(harmony_theta), theta = unname(harmony_theta),
                  reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (res in c(0.3, 0.6, 1.0)) {
  obj <- FindClusters(obj, resolution = res, random.seed = cluster_seed, verbose = FALSE)
}


# Step 5. UMAPs (unintegrated on PCA, integrated on Harmony), seed 42 as rounds 14 / 15
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = umap_seed, reduction.name = "umap_unintegrated", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = umap_seed, verbose = FALSE)


# Step 6. iLISI: for each cell, inverse Simpson index of dataset (and chemistry) among its 30 nearest neighbours
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
say("round 16: %d cells, %d clusters at res 0.6; median iLISI dataset %.2f; T_gd %d", ncol(obj),
    length(unique(obj$RNA_snn_res.0.6)), median(obj$ilisi_dataset), sum(obj$final_label == "T_gd"))


# Step 7. Per-cell table (same columns as cpA14_cells.tsv.gz, plus the round-15 cluster and label)
umap <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
keep_columns <- c("dataset", "sample", "chemistry", "label", "s1_cell_type", "ref_exclude", "r12_cluster", "r12_label", "final_label",
                  "RNA_snn_res.0.3", "RNA_snn_res.0.6", "RNA_snn_res.1", "ilisi_dataset", "ilisi_chemistry",
                  "nCount_RNA", "nFeature_RNA", "nCount_orig", "nFeature_orig", "percent.mt",
                  "r15_cluster_res0.6", "r15_label")
keep_columns <- intersect(keep_columns, names(obj@meta.data))
cells <- data.table(obj@meta.data[, keep_columns], keep.rownames = "cell")
cells$UMAP_1  <- umap[, 1]
cells$UMAP_2  <- umap[, 2]
cells$UMAPu_1 <- umap_unintegrated[, 1]
cells$UMAPu_2 <- umap_unintegrated[, 2]
fwrite(cells, out_file("cells.tsv.gz"), sep = "\t")


# Step 8. How the new clusters (res 0.6) line up with the labels, and with the round-15 clusters
by_cluster <- cells[, .N, by = .(cluster = RNA_snn_res.0.6, final_label)]
by_cluster[, pct := round(100 * N / sum(N), 1), by = cluster]
by_cluster <- by_cluster[order(cluster, -N)]
fwrite(by_cluster, out_file("cluster_by_final_label.tsv"), sep = "\t")
print(by_cluster[, .SD[1], by = cluster][order(pct)])
print(by_cluster[final_label == "T_gd"][order(-N)])

old_vs_new <- cells[, .N, by = .(r15_cluster_res0.6, new_cluster = RNA_snn_res.0.6)]
print(old_vs_new[order(r15_cluster_res0.6, -N)][, .SD[1], by = r15_cluster_res0.6])


# Step 9. Label file for step C, in the same layout as 2026-10-23_round14_labels.tsv.gz
majority <- by_cluster[, .(cluster_type = final_label[1]), by = cluster]
label_file <- cells[, .(cell, dataset, published_label = label, cluster = RNA_snn_res.0.6, final_label)]
label_file$cluster_type <- majority$cluster_type[match(label_file$cluster, majority$cluster)]
setcolorder(label_file, c("cell", "dataset", "published_label", "cluster", "cluster_type", "final_label"))
label_path <- file.path(out_dir, "stepC/inputs/2026-10-26_round16_labels.tsv")
fwrite(label_file, label_path, sep = "\t")
system(sprintf("gzip -f %s", label_path))


saveRDS(obj, file.path(out_dir, "cHL_merged_harmony_v2_r16.rds"))
say("saved round 16")


# ---------------------------------------------------------------------------------------------------------------------
# Part B. Unstable-cell check on round 16 (report only; the saved object is not changed)
# ---------------------------------------------------------------------------------------------------------------------

# Step 10. Recluster at res 0.6 with the new seed (same genes, same recipe)
check_obj <- obj
check_obj@graphs <- list()
check_obj <- RunPCA(check_obj, features = hvg, npcs = 50, seed.use = check_seed, verbose = FALSE)
set.seed(check_seed)
check_obj <- RunHarmony(check_obj, group.by.vars = names(harmony_theta), theta = unname(harmony_theta),
                        reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
check_obj <- FindNeighbors(check_obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
check_obj <- FindClusters(check_obj, resolution = 0.6, random.seed = check_seed, verbose = FALSE)


# Step 11. Cluster type = most common round-16 label in each new cluster; unstable = label != cluster type
check_cells <- data.table(cell = colnames(check_obj), new_cluster = as.character(check_obj$RNA_snn_res.0.6),
                          label = check_obj$final_label, r16_cluster = as.character(obj$RNA_snn_res.0.6),
                          dataset = check_obj$dataset, sample = check_obj$sample, published = check_obj$label,
                          n_umi = check_obj$nCount_RNA)
label_counts <- check_cells[, .N, by = .(new_cluster, label)]
label_counts <- label_counts[order(new_cluster, -N)]
cluster_types <- label_counts[, .(cluster_type = label[1],
                                  type_pct     = round(100 * N[1] / sum(N), 1),
                                  cells        = sum(N),
                                  second       = paste(label[2], N[2])),
                              by = new_cluster]
check_cells <- merge(check_cells, cluster_types[, .(new_cluster, cluster_type)], by = "new_cluster")
check_cells$is_protected <- check_cells$label %in% protected
check_cells$mismatch <- check_cells$label != check_cells$cluster_type
check_cells$unstable <- check_cells$mismatch & !check_cells$is_protected
n_unstable <- sum(check_cells$unstable)
n_protected_mismatch <- sum(check_cells$mismatch & check_cells$is_protected)
say("round-16 check (seed %d): %d cells, %d clusters, %d unstable (%.2f %%); protected cells outside a cluster of their type: %d",
    check_seed, nrow(check_cells), nrow(cluster_types), n_unstable, 100 * n_unstable / nrow(check_cells), n_protected_mismatch)


# Step 12. Summary tables and the flagged cells
by_type <- check_cells[, .(cells = .N, unstable = sum(unstable), protected_mismatch = sum(mismatch & is_protected)), by = label]
by_type$pct_unstable <- round(100 * by_type$unstable / by_type$cells, 2)
by_type <- by_type[order(-unstable)]
moves <- check_cells[unstable == TRUE, .N, by = .(label, cluster_type)]
moves <- moves[order(-N)]
by_dataset <- check_cells[, .(cells = .N, unstable = sum(unstable)), by = dataset]
by_dataset$pct_unstable <- round(100 * by_dataset$unstable / by_dataset$cells, 2)
by_dataset <- by_dataset[order(dataset)]
cluster_types <- cluster_types[order(as.integer(new_cluster))]
print(cluster_types)
print(by_type)
print(moves)
print(by_dataset)

fwrite(data.table(seed = check_seed, cells = nrow(check_cells), clusters = nrow(cluster_types), unstable = n_unstable,
                  protected_mismatch = n_protected_mismatch), check_file("summary.tsv"), sep = "\t")
fwrite(by_type, check_file("by_type.tsv"), sep = "\t")
fwrite(moves, check_file("moves.tsv"), sep = "\t")
fwrite(by_dataset, check_file("by_dataset.tsv"), sep = "\t")
fwrite(cluster_types, check_file("cluster_types.tsv"), sep = "\t")
flagged_r16 <- check_cells[mismatch == TRUE, .(cell, dataset, sample, published, label, cluster_type, new_cluster, r16_cluster,
                                              n_umi, is_protected, unstable)]
fwrite(flagged_r16, check_file("unstable_cells.tsv.gz"), sep = "\t")
say("done")
