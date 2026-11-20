# =============================================================================
# v2_C14_unstable_loop_r13b.R
#
# Purpose (author, 2026-10-23): "Run the unstable loop with 4,000 HVG, protect T_gd."
#   Start: the round-13b object (v2_C13_recluster_after_tgd_relabel.R; 112,364 cells) with the T_gd-sub-UMAP labels.
#   Loop (author pipeline, step 5; same as round 13, v2_C8_round13_loop.R, except 4,000 variable genes):
#     each round: recluster the whole object with a NEW seed (1400 + round number)
#       4,000 variable genes picked per dataset, ranked by how many datasets pick them -> scale -> PCA 50
#       -> Harmony on dataset + chemistry + sample (theta 4 / 4 / 2) -> neighbours on 30 dims -> Louvain res 0.6;
#     cluster type = the majority of its cells' current labels (cell-type level: Treg stays separate from CD4);
#     unstable cell = its label differs from its cluster's type -> removed;
#     protected: T_gd (never removed, keeps its own label; T_gd forms no cluster of its own);
#     kept cells roll their label forward (= unchanged, because a kept unprotected cell already matches its cluster);
#     stop when a round removes fewer than 100 cells (at most 10 rounds).
#   Then a final recluster (seed 1499; res 0.3 / 0.6 / 1.0) + UMAPs (seed 42) + iLISI -> round 14.
#
# Inputs (OAK):  output/cHL_merged_harmony_v2_r13b.rds (final_label = T_gd-sub-UMAP labels)
# Outputs (OAK): output/cHL_merged_harmony_v2_r14.rds
#                output/cpA14_loop_rounds.tsv        (cells in, clusters, removed, % removed per dataset, per round)
#                output/cpA14_removed_cells.tsv.gz   (cell, round, dataset, sample, label, cluster type)
#                output/cpA14_cells.tsv.gz           (per-cell table, same columns as cpA13b_cells.tsv.gz)
#                output/cpA14_cluster_by_final_label.tsv, output/cpA14_hvg_4000.tsv
#                output/stepC/inputs/2026-10-23_round14_labels.tsv.gz (cell, dataset, final_label, cluster)
# Run on Sherlock:  sbatch v2_run.sh v2_C14_unstable_loop_r13b.R
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

n_hvg         <- 4000                                        # author 2026-10-23
harmony_theta <- c(dataset = 4, chemistry = 4, sample = 2)   # as round 13
protected     <- "T_gd"                                      # author: protect T_gd
min_removed   <- 100                                         # stop when a round removes fewer cells than this
max_rounds    <- 10
gene_table <- fread(file.path(work_dir, "ref/gene_table_v47.tsv"))


# Helper: recluster the whole object (variable genes -> PCA -> Harmony -> neighbours -> Louvain).
#   resolutions = 0.6 inside the loop; 0.3 / 0.6 / 1.0 for the final recluster. Returns the object; the chosen
#   variable genes are stored in attr(obj, "hvg").
recluster <- function(obj, seed, resolutions) {
  obj@reductions <- list()
  obj@graphs <- list()
  obj@neighbors <- list()

  # variable genes: n_hvg per dataset, ranked by how many datasets pick them (ties: median rank across datasets)
  hvg_per_dataset <- list()
  for (ds in unique(obj$dataset)) {
    ds_cells  <- colnames(obj)[obj$dataset == ds]
    ds_object <- CreateSeuratObject(LayerData(obj, assay = "RNA", layer = "counts")[, ds_cells])
    ds_object <- FindVariableFeatures(ds_object, nfeatures = n_hvg, verbose = FALSE)
    hvg_per_dataset[[ds]] <- VariableFeatures(ds_object)
  }
  all_hvg <- names(table(unlist(hvg_per_dataset)))
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
  obj <- ScaleData(obj, features = hvg, verbose = FALSE)
  obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = seed, verbose = FALSE)
  set.seed(seed)
  obj <- RunHarmony(obj, group.by.vars = names(harmony_theta), theta = unname(harmony_theta),
                    reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
  for (res in resolutions) {
    obj <- FindClusters(obj, resolution = res, random.seed = seed, verbose = FALSE)
  }
  attr(obj, "hvg") <- hvg
  obj
}


# Step 1. Start from round 13b; every cell's current label = its T_gd-sub-UMAP label
obj <- readRDS(file.path(out_dir, "cHL_merged_harmony_v2_r13b.rds"))
obj$current_label <- obj$final_label
obj$is_protected  <- obj$final_label %in% protected
say("round 13b: %d cells; protected T_gd: %d", ncol(obj), sum(obj$is_protected))


# Step 2. The unstable-cell loop
round_summaries <- list()
removed_tables  <- list()
round <- 0
repeat {
  round <- round + 1
  obj <- recluster(obj, seed = 1400 + round, resolutions = 0.6)

  cells <- data.table(cell = colnames(obj), cluster = as.character(obj$RNA_snn_res.0.6), label = obj$current_label,
                      is_protected = obj$is_protected, dataset = obj$dataset, sample = obj$sample)

  # cluster type = majority current label of the cluster's cells
  label_counts <- cells[, .N, by = .(cluster, label)]
  label_counts <- label_counts[order(cluster, -N)]
  cluster_types <- label_counts[, .(cluster_type = label[1]), by = cluster]
  cells <- merge(cells, cluster_types, by = "cluster")

  cells$unstable <- cells$label != cells$cluster_type & !cells$is_protected
  n_removed <- sum(cells$unstable)
  n_protected_mismatch <- sum(cells$label != cells$cluster_type & cells$is_protected)
  say("round %d: %d cells, %d clusters, %d unstable removed (%d protected T_gd outside a T_gd cluster, kept)",
      round, nrow(cells), nrow(cluster_types), n_removed, n_protected_mismatch)
  swaps <- cells[unstable == TRUE, .N, by = .(label, cluster_type)]
  print(head(swaps[order(-N)], 12))

  # per-round summary, with % removed per dataset
  summary_row <- data.table(round = round, cells_in = nrow(cells), clusters = nrow(cluster_types),
                            removed = n_removed, protected_mismatch = n_protected_mismatch)
  for (ds in sort(unique(cells$dataset))) {
    summary_row[[paste0("pct_removed_", ds)]] <- round(100 * mean(cells$unstable[cells$dataset == ds]), 2)
  }
  round_summaries[[round]] <- summary_row
  removed_tables[[round]] <- cells[unstable == TRUE, .(cell, round = round, dataset, sample, label, cluster_type)]

  # remove unstable cells; kept unprotected cells take their cluster's type (rolling; no change for them)
  kept <- cells[unstable == FALSE]
  obj <- subset(obj, cells = kept$cell)
  kept <- kept[match(colnames(obj), kept$cell)]
  roll <- !kept$is_protected
  obj$current_label[roll] <- kept$cluster_type[roll]

  if (n_removed < min_removed || round >= max_rounds) break
}
if (n_removed >= min_removed) cat("WARNING: stopped at", max_rounds, "rounds before a round removed fewer than", min_removed, "cells\n")

round_summaries <- rbindlist(round_summaries, fill = TRUE)
print(round_summaries)
fwrite(round_summaries, file.path(out_dir, "cpA14_loop_rounds.tsv"), sep = "\t")
fwrite(rbindlist(removed_tables), file.path(out_dir, "cpA14_removed_cells.tsv.gz"), sep = "\t")


# Step 3. Final recluster (seed 1499) at res 0.3 / 0.6 / 1.0, UMAPs (seed 42)
obj$final_label <- obj$current_label
obj <- recluster(obj, seed = 1499, resolutions = c(0.3, 0.6, 1.0))
hvg <- attr(obj, "hvg")
fwrite(data.table(ensg = hvg, symbol = gene_table$symbol[match(hvg, gene_table$ensg)]), file.path(out_dir, "cpA14_hvg_4000.tsv"), sep = "\t")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 42, reduction.name = "umap_unintegrated", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)


# Step 4. iLISI: inverse Simpson index of dataset (and chemistry) among each cell's 30 nearest neighbours
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
say("round 14: %d cells, %d clusters at res 0.6; median iLISI dataset %.2f; T_gd %d", ncol(obj),
    length(unique(obj$RNA_snn_res.0.6)), median(obj$ilisi_dataset), sum(obj$final_label == "T_gd"))


# Step 5. Exports: per-cell table, cluster x label table, label file for step C, object
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
fwrite(cells, file.path(out_dir, "cpA14_cells.tsv.gz"), sep = "\t")

by_cluster <- cells[, .N, by = .(cluster = RNA_snn_res.0.6, final_label)]
by_cluster[, pct := round(100 * N / sum(N), 1), by = cluster]
by_cluster <- by_cluster[order(cluster, -N)]
fwrite(by_cluster, file.path(out_dir, "cpA14_cluster_by_final_label.tsv"), sep = "\t")
print(by_cluster[, .SD[1], by = cluster][order(pct)])

# label file in the same layout as 2026-10-21_round13_labels.tsv.gz (step C needs final_label, cluster, cluster_type)
majority <- by_cluster[, .(cluster_type = final_label[1]), by = cluster]
label_file <- cells[, .(cell, dataset, published_label = label, cluster = RNA_snn_res.0.6, final_label)]
label_file$cluster_type <- majority$cluster_type[match(label_file$cluster, majority$cluster)]
setcolorder(label_file, c("cell", "dataset", "published_label", "cluster", "cluster_type", "final_label"))
fwrite(label_file, file.path(out_dir, "stepC/inputs/2026-10-23_round14_labels.tsv"), sep = "\t")
system(sprintf("gzip -f %s", file.path(out_dir, "stepC/inputs/2026-10-23_round14_labels.tsv")))

saveRDS(obj, file.path(out_dir, "cHL_merged_harmony_v2_r14.rds"))
say("saved round 14")
