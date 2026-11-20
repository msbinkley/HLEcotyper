# =============================================================================
# v2_C8_round13_loop.R — round 13 of the merged cHL object (runs on Sherlock).
#
# PURPOSE (author 2026-10-20)
#   1. Remove the 275 cells of the cluster-8 (T_MKI67) sub-UMAP subclusters 10 and 12 (low UMI; v2_C7 round 1, res 0.8)
#      from r12.
#   2. Remove unstable cells until a round removes < 100 cells (author pipeline, rolling annotation, cell-type level):
#      - every cell starts with its round-12 label (r12 cluster annotation + sub-UMAP T_gd; 2026-10-20_round12_labels.tsv.gz),
#        with T_MKI67 renamed "Cycling" (r12 cluster 8 = cycling T + cycling GC B; author 2026-10-20);
#      - each round: recluster with a new seed (per-dataset HVG 3,000 -> PCA 50 -> Harmony dataset + chemistry + sample,
#        theta 4/4/2 -> Louvain res 0.6); cluster type = majority of its cells' current labels; unstable = label !=
#        cluster type -> removed; kept cells roll their label forward (= unchanged for unprotected cells). Max 10 rounds.
#      - Protected: T_gd only (never forms its own cluster) - never removed, keeps its label. HRS / ILC form their own
#        clusters (r12 17 / 15) and are not protected (author 2026-10-20).
#   3. Final full recluster + UMAP -> round 13.
#
# INPUT  : output/cHL_merged_harmony_v2_r12.rds, output/stepC/inputs/2026-10-20_round12_labels.tsv.gz,
#          output/subumap/r12_cl8_r2_removed_cells.tsv.gz, ref/gene_table_v47.tsv
# OUTPUT : output/cHL_merged_harmony_v2_r13.rds, output/cpA13_* (cells, rounds, removed cells, HVG, dot plots, cluster check)
# RUN    : on Sherlock (sbatch v2_run.sh v2_C8_round13_loop.R)
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
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))

# Print a time-stamped log line
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

MIN_DROP   <- 100L                                          # stop when a round removes fewer than 100 cells
MAX_ROUNDS <- 10L                                           # ... or after 10 rounds
THETA      <- c(dataset = 4, chemistry = 4, sample = 2)     # Harmony theta per batch variable
PROTECT    <- "T_gd"                                        # protected type: never removed, keeps its label


# Helper: recluster the whole object with the given seed.
#   per-dataset HVG 3,000 -> scale -> PCA 50 -> Harmony (dataset + chemistry + sample, theta 4/4/2) -> neighbours on
#   Harmony dims 1-30 -> Louvain res 0.6 (final = TRUE: res 0.3, 0.6 and 1.0). The chosen HVGs are attached as attr "hvg".
recluster <- function(obj, seed, final = FALSE) {
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
  if (final) {
    resolutions <- c(0.3, 0.6, 1.0)
  } else {
    resolutions <- 0.6
  }
  for (r in resolutions) {
    obj <- FindClusters(obj, resolution = r, random.seed = seed, verbose = FALSE)
  }
  attr(obj, "hvg") <- hvg
  obj
}


## Step 1. r12 minus the 275 cluster-8 low-UMI cells; round-12 labels ---------------------------------------------------------
obj   <- readRDS(file.path(OUT, "cHL_merged_harmony_v2_r12.rds"))
lab   <- fread(cmd = sprintf("gzip -dc %s", file.path(OUT, "stepC/inputs/2026-10-20_round12_labels.tsv.gz")))
drop8 <- fread(cmd = sprintf("gzip -dc %s", file.path(OUT, "subumap/r12_cl8_r2_removed_cells.tsv.gz")))
stopifnot(nrow(drop8) == 275, all(drop8$cell %in% colnames(obj)), setequal(lab$cell, colnames(obj)))
obj$r12_cluster <- as.character(obj$RNA_snn_res.0.6)
obj$r12_label   <- lab$r12_label[match(colnames(obj), lab$cell)]
## author 2026-10-20: r12 cluster 8 is cycling T + cycling GC B -> its type is "Cycling" from now on (was T_MKI67). All 4,198
## T_MKI67-labelled r12 cells are cluster 8; in the loop they are stable while they stay in a Cycling-majority cluster.
stopifnot(all(obj$r12_cluster[obj$r12_label == "T_MKI67"] == "8"))
obj$r12_label[obj$r12_label == "T_MKI67"] <- "Cycling"
obj <- subset(obj, cells = setdiff(colnames(obj), drop8$cell))
obj$cell_type <- obj$r12_label                       # current label (rolled forward each round)
obj$loop_keep <- obj$r12_label %in% PROTECT          # protected cells (T_gd)
log_("r12 %d cells -> removed %d (cluster-8 sub-UMAP 10 / 12) -> %d; protected T_gd %d",
     length(lab$cell), nrow(drop8), ncol(obj), sum(obj$loop_keep))


## Step 2. Unstable-cell loop ------------------------------------------------------------------------------------------------
rounds  <- list()            # one summary row per round
removed <- list()            # removed cells of each round
n <- 0L                      # round number
repeat {
  n <- n + 1L
  obj <- recluster(obj, seed = 1300L + n)            # new seed every round: 1301, 1302, ...

  # one row per cell: new cluster, current label, protected flag
  loop_table <- data.table(cell = colnames(obj), cl = as.character(obj$RNA_snn_res.0.6), type = obj$cell_type,
                           keep = obj$loop_keep, dataset = obj$dataset, sample = obj$sample, label = obj$label)
  # cluster type = most common current label in the new cluster (and its %)
  type_counts <- loop_table[, .N, by = .(cl, type)]
  type_counts <- type_counts[order(cl, -N)]
  ct <- type_counts[, .(cluster_type = type[1], type_pct = round(100 * N[1] / sum(N), 1)), by = cl]
  loop_table <- merge(loop_table, ct, by = "cl")
  # unstable = label differs from the cluster type, unless protected
  loop_table[, unstable := type != cluster_type & !keep]
  n_drop <- sum(loop_table$unstable)
  n_protected_mismatch <- loop_table[type != cluster_type & keep, .N]
  log_("loop round %d: %d cells, %d clusters, %d unstable (+%d protected T_gd would have been unstable)",
       n, nrow(loop_table), nrow(ct), n_drop, n_protected_mismatch)
  moves <- loop_table[unstable == TRUE, .N, by = .(type, cluster_type)]
  print(head(moves[order(-N)], 12))

  # round summary: counts + % of each dataset's cells removed (columns pct_drop_<dataset>)
  pct <- loop_table[, .(pct = round(100 * mean(unstable), 2)), by = dataset]
  pct <- pct[order(dataset)]
  round_info <- list(round = n, cells_in = nrow(loop_table), clusters = nrow(ct), dropped = n_drop,
                     protected_mismatch = n_protected_mismatch)
  pct_info <- as.list(setNames(pct$pct, paste0("pct_drop_", pct$dataset)))
  rounds[[n]] <- as.data.table(c(round_info, pct_info))
  removed[[n]] <- loop_table[unstable == TRUE, .(cell, round = n, dataset, sample, label, type, cluster_type)]

  # remove unstable cells; kept unprotected cells take their cluster's type (rolling; no-op for kept cells)
  obj <- subset(obj, cells = loop_table[unstable == FALSE, cell])
  kept_rows <- loop_table[match(colnames(obj), cell)]
  roll <- !kept_rows$keep
  obj$cell_type[roll] <- kept_rows$cluster_type[roll]

  if (n_drop < MIN_DROP || n >= MAX_ROUNDS) break
}
if (n_drop >= MIN_DROP) cat("WARNING: stopped at", MAX_ROUNDS, "rounds before a round removed fewer than", MIN_DROP, "cells\n")
round_summary <- rbindlist(rounds, fill = TRUE)
print(round_summary)
fwrite(round_summary, file.path(OUT, "cpA13_loop_rounds.tsv"), sep = "\t")

# all removed cells: round 0 = the 275 cluster-8 cells, then each loop round
removed_cluster8 <- drop8[, .(cell, round = 0L, reason = "cluster-8 sub-UMAP subcluster 10/12 (low UMI)")]
removed_loop <- rbindlist(removed)
removed_loop <- removed_loop[, .(cell, round, reason = paste0("unstable: ", type, " in ", cluster_type, " cluster"))]
fwrite(rbind(removed_cluster8, removed_loop), file.path(OUT, "cpA13_removed_cells.tsv.gz"), sep = "\t")


## Step 3. Final recluster + UMAP + exports ----------------------------------------------------------------------------------
obj$r13_label <- obj$cell_type
obj <- recluster(obj, seed = 1399L, final = TRUE)
hvg <- attr(obj, "hvg")
fwrite(data.table(ensg = hvg, symbol = v47$symbol[match(hvg, v47$ensg)]), file.path(OUT, "cpA13_hvg_3000.tsv"), sep = "\t")
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 42, reduction.name = "umap_unintegrated", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)

# Mixing score (iLISI): per cell, inverse Simpson index of dataset / chemistry among its 30 nearest neighbours
neighbours <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, k.param = 30, return.neighbor = TRUE, verbose = FALSE)@neighbors[[1]]
neighbour_index <- Indices(neighbours)
inverse_simpson <- function(lab) {
  apply(neighbour_index, 1, function(r) {
    p <- table(lab[r]) / length(r)                               # share of each group among the neighbours
    1 / sum(p^2)
  })
}
obj$ilisi_dataset   <- inverse_simpson(obj$dataset)
obj$ilisi_chemistry <- inverse_simpson(obj$chemistry)
log_("round 13: %d cells, %d clusters (res 0.6); median iLISI dataset %.2f (S1 %.2f); T_gd %d", ncol(obj), length(unique(obj$RNA_snn_res.0.6)),
     median(obj$ilisi_dataset), median(obj$ilisi_dataset[obj$dataset == "s1"]), sum(obj$r13_label == "T_gd"))

# Per-cell export with UMAP coordinates
umap_integrated   <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
cl <- grep("^RNA_snn_res", names(obj@meta.data), value = TRUE)
keep_cols <- intersect(c("dataset", "sample", "chemistry", "label", "s1_cell_type", "ref_exclude", "r12_cluster", "r12_label", "r13_label", cl,
                         "ilisi_dataset", "ilisi_chemistry", "nCount_RNA", "nFeature_RNA", "nCount_orig", "nFeature_orig", "percent.mt"), names(obj@meta.data))
cells <- data.table(obj@meta.data[, keep_cols], keep.rownames = "cell")
cells[, `:=`(UMAP_1  = umap_integrated[, 1],   UMAP_2  = umap_integrated[, 2],
             UMAPu_1 = umap_unintegrated[, 1], UMAPu_2 = umap_unintegrated[, 2])]
setnames(cells, "r13_label", "final_label")                                   # v2_cpA_plots.R draws "final_label" as its third panel
fwrite(cells, file.path(OUT, "cpA13_cells.tsv.gz"), sep = "\t")

# Each res-0.6 cluster: cells per final label and their % of the cluster
chk <- cells[, .N, by = .(cluster = RNA_snn_res.0.6, final_label)]
chk[, pct := round(100 * N / sum(N), 1), by = cluster]
chk <- chk[order(cluster, -N)]
fwrite(chk, file.path(OUT, "cpA13_cluster_by_final_label.tsv"), sep = "\t")
print(chk[, .SD[1], by = cluster][order(pct)])                   # main label of each cluster, least pure first

# Dot-plot tables: per group, mean log-expression and % cells expressing each marker gene
markers <- c("PTPRC","CD3E","CD4","CD8A","FOXP3","IL2RA","MKI67","TOP2A","NKG7","GNLY","KLRD1","TRDC","TRGC1","KIT","IL7R","MS4A1","CD19","TCL1A",
             "IGHD","CD27","BCL6","AICDA","RGS13","MZB1","JCHAIN","SDC1","CD68","CD163","C1QA","LYZ","CD14","VCAN","FCGR3A","CLEC9A","CD1C","LILRA4",
             "IL3RA","COL1A1","COL1A2","PDGFRA","PECAM1","VWF","TNFRSF8","IRF4","CD274","FOS","JUN","HSPA1A","EGR1")
mg <- v47[symbol %in% markers & ensg %in% rownames(obj)]           # marker genes present in the object (v47 table order)
dat <- LayerData(obj, assay = "RNA", layer = "data")[mg$ensg, ]
dot <- function(group) {
  cell_index_by_group <- split(seq_len(ncol(obj)), group)
  per_group <- vector("list", length(cell_index_by_group))
  names(per_group) <- names(cell_index_by_group)
  for (k in seq_along(cell_index_by_group)) {
    i <- cell_index_by_group[[k]]
    per_group[[k]] <- data.table(gene        = mg$symbol,
                                 mean_logexp = Matrix::rowMeans(dat[, i, drop = FALSE]),
                                 pct_expr    = 100 * Matrix::rowMeans(dat[, i, drop = FALSE] > 0))
  }
  rbindlist(per_group, idcol = "group")
}
fwrite(dot(obj$RNA_snn_res.0.6), file.path(OUT, "cpA13_dotplot_by_cluster_res0.6.tsv"), sep = "\t")
fwrite(dot(obj$r13_label), file.path(OUT, "cpA13_dotplot_by_final_label.tsv"), sep = "\t")
# by published label; S1 cells (no published label) by "S1:<ex vivo type>"
fwrite(dot(ifelse(is.na(obj$label), paste0("S1:", obj$s1_cell_type), obj$label)), file.path(OUT, "cpA13_dotplot_by_label.tsv"), sep = "\t")
saveRDS(obj, file.path(OUT, "cHL_merged_harmony_v2_r13.rds"))
log_("saved round 13")
