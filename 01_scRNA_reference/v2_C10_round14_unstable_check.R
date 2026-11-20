# =============================================================================
# v2_C10_round14_unstable_check.R — one more unstable-cell round on the final round-13 object (runs on Sherlock).
#
# PURPOSE (author 2026-10-23)
#   "for the cHL sc/sn RNA-seq reference (merged), perform another round of unstable cell identification, report how many
#   cells are removed".
#   Same round as v2_C8_round13_loop.R (cell-type level, rolling labels):
#     recluster with a NEW seed (1401; r13 used 1301-1305 and 1399) -> per-dataset HVG 3,000 -> PCA 50 -> Harmony dataset +
#     chemistry + sample, theta 4/4/2 -> Louvain res 0.6; cluster type = majority of its cells' round-13 final labels;
#     unstable = label != cluster type. Protected: T_gd (as round 13).
#   Every cell is included whatever its source dataset (incl. the 157 10x generic-B cells, with their final labels).
#   Report only: the object is NOT changed or saved; removal is the author's decision.
#
# INPUT  : output/cHL_merged_harmony_v2_r13.rds, output/stepC/inputs/2026-10-21_round13_labels.tsv.gz
# OUTPUT : output/cpA14check_{summary,by_type,moves,by_dataset,cluster_types}.tsv, output/cpA14check_unstable_cells.tsv.gz
# RUN    : on Sherlock (sbatch v2_run.sh v2_C10_round14_unstable_check.R)
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
})
W   <- "${PROJECT_ROOT}/cHL_v2"      # project root on OAK
OUT <- file.path(W, "output")

# Print a time-stamped log line
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

SEED    <- 1401L                                         # new seed (round 13 used 1301-1305 and 1399)
THETA   <- c(dataset = 4, chemistry = 4, sample = 2)     # Harmony theta per batch variable
PROTECT <- "T_gd"                                        # protected type: never counted as unstable


# Helper: recluster the whole object = v2_C8_round13_loop.R recluster() with final = FALSE.
#   per-dataset HVG 3,000 -> scale -> PCA 50 -> Harmony (dataset + chemistry + sample, theta 4/4/2) -> neighbours on
#   Harmony dims 1-30 -> Louvain res 0.6
recluster <- function(obj, seed) {
  obj@reductions <- list()
  obj@graphs     <- list()
  obj@neighbors  <- list()

  # 3,000 variable genes per dataset
  cells_by_dataset <- split(colnames(obj), obj$dataset)
  hv <- list()
  for (ds in names(cells_by_dataset)) {
    ds_obj <- CreateSeuratObject(LayerData(obj, assay = "RNA", layer = "counts")[, cells_by_dataset[[ds]]])
    ds_obj <- FindVariableFeatures(ds_obj, nfeatures = 3000, verbose = FALSE)
    hv[[ds]] <- VariableFeatures(ds_obj)
  }
  rm(ds_obj)

  # rank genes by how many datasets pick them, then by their median rank within those datasets; keep the top 3,000
  rk <- table(unlist(hv))
  rank_matrix <- matrix(NA_integer_, nrow = length(rk), ncol = length(hv))   # gene x dataset: rank in that dataset's list
  for (d in seq_along(hv)) {
    rank_matrix[, d] <- match(names(rk), hv[[d]])
  }
  med <- apply(rank_matrix, 1, median, na.rm = TRUE)
  hvg <- names(rk)[order(-rk, med)][1:3000]
  DefaultAssay(obj) <- "RNA"
  VariableFeatures(obj) <- hvg

  obj <- ScaleData(obj, features = hvg, verbose = FALSE)
  obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = seed, verbose = FALSE)
  set.seed(seed)                                                # Harmony uses R's random numbers
  obj <- RunHarmony(obj, group.by.vars = names(THETA), theta = unname(THETA), reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
  FindClusters(obj, resolution = 0.6, random.seed = seed, verbose = FALSE)
}


## Step 1. Round-13 object and its final labels (must agree) -------------------------------------------------------------------
obj <- readRDS(file.path(OUT, "cHL_merged_harmony_v2_r13.rds"))
lab <- fread(cmd = sprintf("gzip -dc %s", file.path(OUT, "stepC/inputs/2026-10-21_round13_labels.tsv.gz")))
stopifnot(setequal(lab$cell, colnames(obj)),
          all(obj$r13_label == lab$final_label[match(colnames(obj), lab$cell)]))
log_("round 13: %d cells; T_gd %d (protected)", ncol(obj), sum(obj$r13_label %in% PROTECT))


## Step 2. Recluster with the new seed and find unstable cells ------------------------------------------------------------------
obj <- recluster(obj, SEED)
cell_table <- data.table(cell = colnames(obj), cl = as.character(obj$RNA_snn_res.0.6), type = obj$r13_label,
                         keep = obj$r13_label %in% PROTECT, dataset = obj$dataset, sample = obj$sample, published = obj$label)

# cluster type = most common round-13 label in the new cluster; also its %, cluster size, and the runner-up label
type_counts <- cell_table[, .N, by = .(cl, type)]
type_counts <- type_counts[order(cl, -N)]
ct <- type_counts[, .(cluster_type = type[1],
                      type_pct     = round(100 * N[1] / sum(N), 1),
                      cells        = sum(N),
                      second       = paste(type[2], N[2])),
                  by = cl]
cell_table <- merge(cell_table, ct[, .(cl, cluster_type)], by = "cl")
cell_table[, unstable := type != cluster_type & !keep]          # unstable = label differs from cluster type, unless protected
n_drop <- cell_table[unstable == TRUE, .N]
n_protected_mismatch <- cell_table[type != cluster_type & keep, .N]
log_("round 14 check (seed %d): %d cells, %d clusters, %d unstable (%.2f %%); protected T_gd outside a T_gd-majority cluster: %d",
     SEED, nrow(cell_table), nrow(ct), n_drop, 100 * n_drop / nrow(cell_table), n_protected_mismatch)


## Step 3. Summary tables ---------------------------------------------------------------------------------------------------------
# per label: cells, unstable cells (0 if none), % unstable; most unstable first
cells_per_type    <- cell_table[, .(cells = .N), by = type]
unstable_per_type <- cell_table[unstable == TRUE, .(unstable = .N), by = type]
by_type <- merge(cells_per_type, unstable_per_type, by = "type", all.x = TRUE)
by_type[is.na(unstable), unstable := 0L]
by_type[, pct := round(100 * unstable / cells, 2)]
by_type <- by_type[order(-unstable)]

# label -> cluster type of the unstable cells
moves <- cell_table[unstable == TRUE, .N, by = .(type, cluster_type)]
moves <- moves[order(-N)]

# per dataset: cells, unstable, % unstable
by_ds <- cell_table[, .(cells = .N, unstable = sum(unstable), pct = round(100 * mean(unstable), 2)), by = dataset]
by_ds <- by_ds[order(dataset)]

cluster_types_sorted <- ct[order(as.integer(cl))]
print(cluster_types_sorted)
print(by_type)
print(moves)
print(by_ds)


## Step 4. Write the report files (the object is not modified) ---------------------------------------------------------------------
fwrite(data.table(seed = SEED, cells = nrow(cell_table), clusters = nrow(ct), unstable = n_drop, protected_mismatch = n_protected_mismatch),
       file.path(OUT, "cpA14check_summary.tsv"), sep = "\t")
fwrite(by_type, file.path(OUT, "cpA14check_by_type.tsv"), sep = "\t")
fwrite(moves, file.path(OUT, "cpA14check_moves.tsv"), sep = "\t")
fwrite(by_ds, file.path(OUT, "cpA14check_by_dataset.tsv"), sep = "\t")
fwrite(cluster_types_sorted, file.path(OUT, "cpA14check_cluster_types.tsv"), sep = "\t")
fwrite(cell_table[unstable == TRUE, .(cell, dataset, sample, published, type, cluster_type, cl)],
       file.path(OUT, "cpA14check_unstable_cells.tsv.gz"), sep = "\t")
log_("done (object not modified)")
