# =============================================================================
# v2_C21_round15_unstable_check.R — one unstable-cell round on round 15 (report only; runs on Sherlock).
#
# Purpose (author, 2026-10-26):
#   "Check if there are unstable cells in r15 that should be removed"
#
# One round of the round-14 unstable loop (v2_C14_unstable_loop_r13b.R, step 2), cell-type level:
#   recluster the whole object with a NEW seed (1501; round 14 used 1401-1405 and 1499, round 15 used 1499)
#   -> 4,000 variable genes per dataset, ranked by how many datasets pick them -> scale -> PCA 50
#   -> Harmony on dataset + chemistry + sample (theta 4 / 4 / 2) -> neighbours on 30 Harmony dims -> Louvain res 0.6
#   cluster type = the most common round-15 label in the new cluster;
#   unstable = label != cluster type.
# Protected (never counted as unstable; reported separately): T_gd and HRS (author rule: rare types are protected by
#   their own label).
# Report only: the object is NOT changed or saved; removal is the author's decision.
#
# Input (OAK):  output/cHL_merged_harmony_v2_r15.rds
# Outputs (OAK): output/2026-10-26_cpA15check_{summary,by_type,moves,by_dataset,cluster_types}.tsv
#                output/2026-10-26_cpA15check_unstable_cells.tsv.gz  (unstable cells + protected mismatches, flagged)
# Run on Sherlock:  sbatch v2_run.sh v2_C21_round15_unstable_check.R
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
out_file <- function(name) file.path(out_dir, paste0("2026-10-26_cpA15check_", name))
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

check_seed    <- 1501                                        # new seed (not used by rounds 14 or 15)
n_hvg         <- 4000                                        # as rounds 14 and 15
harmony_theta <- c(dataset = 4, chemistry = 4, sample = 2)   # as rounds 14 and 15
protected     <- c("T_gd", "HRS")                            # author rule: rare types are never removed


# Step 1. Read round 15
obj <- readRDS(file.path(out_dir, "cHL_merged_harmony_v2_r15.rds"))
say("round 15: %d cells; protected: T_gd %d, HRS %d", ncol(obj), sum(obj$final_label == "T_gd"), sum(obj$final_label == "HRS"))
obj$r15_cluster <- as.character(obj$RNA_snn_res.0.6)            # keep the round-15 clusters (FindClusters overwrites res 0.6)


# Step 2. Variable genes: 4,000 per dataset, then rank genes by how many datasets pick them
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


# Step 3. Scale, PCA 50, Harmony, neighbours, Louvain res 0.6 with the new seed
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = check_seed, verbose = FALSE)
set.seed(check_seed)
obj <- RunHarmony(obj, group.by.vars = names(harmony_theta), theta = unname(harmony_theta),
                  reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.6, random.seed = check_seed, verbose = FALSE)


# Step 4. Cluster type = most common round-15 label in each new cluster; unstable = label != cluster type
cells <- data.table(cell = colnames(obj), new_cluster = as.character(obj$RNA_snn_res.0.6), label = obj$final_label,
                    r15_cluster = obj$r15_cluster, dataset = obj$dataset, sample = obj$sample,
                    published = obj$label, n_umi = obj$nCount_RNA)
label_counts <- cells[, .N, by = .(new_cluster, label)]
label_counts <- label_counts[order(new_cluster, -N)]
cluster_types <- label_counts[, .(cluster_type = label[1],
                                  type_pct     = round(100 * N[1] / sum(N), 1),
                                  cells        = sum(N),
                                  second       = paste(label[2], N[2])),
                              by = new_cluster]
cells <- merge(cells, cluster_types[, .(new_cluster, cluster_type)], by = "new_cluster")
cells$is_protected <- cells$label %in% protected
cells$mismatch <- cells$label != cells$cluster_type
cells$unstable <- cells$mismatch & !cells$is_protected
n_unstable <- sum(cells$unstable)
n_protected_mismatch <- sum(cells$mismatch & cells$is_protected)
say("check (seed %d): %d cells, %d clusters, %d unstable (%.2f %%); protected cells outside a cluster of their type: %d",
    check_seed, nrow(cells), nrow(cluster_types), n_unstable, 100 * n_unstable / nrow(cells), n_protected_mismatch)


# Step 5. Summary tables
# 5a. per label: cells, unstable cells, %; protected labels show their mismatches in a separate column
by_type <- cells[, .(cells = .N, unstable = sum(unstable), protected_mismatch = sum(mismatch & is_protected)), by = label]
by_type$pct_unstable <- round(100 * by_type$unstable / by_type$cells, 2)
by_type <- by_type[order(-unstable)]

# 5b. label -> cluster type of the unstable cells
moves <- cells[unstable == TRUE, .N, by = .(label, cluster_type)]
moves <- moves[order(-N)]

# 5c. per dataset
by_dataset <- cells[, .(cells = .N, unstable = sum(unstable)), by = dataset]
by_dataset$pct_unstable <- round(100 * by_dataset$unstable / by_dataset$cells, 2)
by_dataset <- by_dataset[order(dataset)]

cluster_types <- cluster_types[order(as.integer(new_cluster))]
print(cluster_types)
print(by_type)
print(moves)
print(by_dataset)


# Step 6. Write the report (the object is not modified or saved)
fwrite(data.table(seed = check_seed, cells = nrow(cells), clusters = nrow(cluster_types), unstable = n_unstable,
                  protected_mismatch = n_protected_mismatch), out_file("summary.tsv"), sep = "\t")
fwrite(by_type, out_file("by_type.tsv"), sep = "\t")
fwrite(moves, out_file("moves.tsv"), sep = "\t")
fwrite(by_dataset, out_file("by_dataset.tsv"), sep = "\t")
fwrite(cluster_types, out_file("cluster_types.tsv"), sep = "\t")
flagged <- cells[mismatch == TRUE, .(cell, dataset, sample, published, label, cluster_type, new_cluster, r15_cluster,
                                     n_umi, is_protected, unstable)]
fwrite(flagged, out_file("unstable_cells.tsv.gz"), sep = "\t")
say("done (object not modified)")
