# =============================================================================
# v2_C3b_resolution_scan.R — r10: which clustering resolution gives T_gd / ILC their own cluster? (Sherlock)
#
# Re-clusters the whole r10 object at resolutions 1.0, 1.5, 2.0, 2.5, 3.0 and 4.0 (same Harmony graph, seed 0).
# For each resolution and each rare type (T_gd, ILC) it finds the cluster with the highest share of that type and reports:
#   best_cluster, cluster_cells, type_cells, purity_pct (% of the cluster that is the type),
#   captures_pct_of_type (% of all cells of the type that fall in this cluster),
#   cluster_majority (most common published label in the cluster), own_cluster (TRUE when that majority is the type).
# Input : output/cHL_merged_harmony_v2_r10.rds
# Output: output/r10_resolution_scan_Tgd_ILC.tsv (also printed)
# Run   : Rscript v2_C3b_resolution_scan.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
})
W <- "${PROJECT_ROOT}/cHL_v2"             # project folder on Sherlock

# Step 1. Load r10; build the neighbour graph only if the object has none
obj <- readRDS(file.path(W, "output/cHL_merged_harmony_v2_r10.rds"))
if (!length(obj@graphs)) obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
cat("graphs:", names(obj@graphs), "\n")

# Step 2. For each resolution: cluster, then find the best cluster for each rare type
labels <- obj$label
results <- list()
for (resolution in c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0)) {
  clusters <- as.character(FindClusters(obj, resolution = resolution, random.seed = 0, verbose = FALSE)$seurat_clusters)
  cell_table <- data.table(cl = clusters, lab = labels)

  for (rare_type in c("T_gd", "ILC")) {
    # cells per cluster (n) and cells of the type per cluster (k); only clusters that contain the type
    per_cluster <- cell_table[, .(n = .N, k = sum(lab %in% rare_type)), by = cl]
    per_cluster <- per_cluster[k > 0]
    per_cluster[, purity := round(100 * k / n, 1)]             # % of the cluster that is the type
    per_cluster[, capture := round(100 * k / sum(k), 1)]       # % of the type's cells in this cluster
    per_cluster <- per_cluster[order(-purity)]
    best <- per_cluster[1]                                     # the purest cluster

    # most common published label in that cluster (unlabelled cells ignored)
    best_cluster_labels <- cell_table[cl == best$cl & !is.na(lab), lab]
    majority <- names(which.max(table(best_cluster_labels)))

    results[[length(results) + 1]] <- data.table(res = resolution, clusters = uniqueN(clusters), type = rare_type,
                                                 best_cluster = best$cl, cluster_cells = best$n,
                                                 type_cells = best$k, purity_pct = best$purity,
                                                 captures_pct_of_type = best$capture,
                                                 cluster_majority = majority, own_cluster = majority == rare_type)
  }
  cat("res", resolution, "done:", uniqueN(clusters), "clusters\n")
}

# Step 3. Save and print
out <- rbindlist(results)
print(out)
fwrite(out, file.path(W, "output/r10_resolution_scan_Tgd_ILC.tsv"), sep = "\t")
