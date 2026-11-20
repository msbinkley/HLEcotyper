# =============================================================================
# v2_C3_unstable_loop.R — unstable-cell loop after the r3 filtering (Sherlock).
#
# Author decision 2026-10-14: remove unstable cells between filterings, until a round removes fewer than 200 cells.
# Rule = author's ex vivo 04_unstable_loop.R, compared at CELL-TYPE level (author 2026-10-14; a first run compared lineages,
# removed 46 cells and is kept in output/superseded_lineage_loop/):
#   - every cell gets the cell type of its r3 cluster (that cluster's majority published label);
#   - recluster with a new seed -> type of each new cluster = the cell type most of its cells carry -> a cell is unstable
#     when its cell type differs from its cluster's type -> drop the unstable cells -> repeat until a round drops < MIN_DROP cells.
# Always kept (author 2026-10-13, rare types): cells of the HRS cluster and every cell published as HRS or T_gd
#   (and ILC, author 2026-10-14; see Settings).
# Recluster = rounds 1-3 recipe (per-dataset HVG -> PCA 50 -> Harmony dataset + chemistry -> res 0.6) with seed = round number.
#
# Settings (environment variables):
#   MIN_DROP     stop when a round drops fewer cells than this (default 200)
#                author 2026-10-14: first 200 (r3 -> r7), then "until < 100" (continued from r7 with START_ROUND=7 MIN_DROP=100)
#   START_ROUND  round to start from (default 3 = the r3 object; otherwise cHL_merged_harmony_v2_r<START_ROUND>.rds)
#   ROLLING      "1" = after each round every kept cell is re-annotated with its new cluster's type (default "0")
#   at most MAX_ROUNDS = 10 rounds per run.
#
# Input : output/cHL_merged_harmony_v2_r3.rds (or the START_ROUND object)
# Output: output/loop/cell_ids_rN.rds per round, output/cpA_loop_rounds.tsv, output/cpA_loop_removed_cells.tsv.gz,
#         final round: output/cHL_merged_harmony_v2_<rN>.rds + output/cpA<N>_* tables (same layout as cpA3_, for v2_cpA_plots.R)
#         (file names get "_from_r<START_ROUND>" when START_ROUND != 3 and "_rolling" / "rolling" with ROLLING=1,
#          so each run keeps its own tables)
# Run   : [MIN_DROP=100] [START_ROUND=7] [ROLLING=1] Rscript v2_C3_unstable_loop.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
})
W <- "${PROJECT_ROOT}/cHL_v2"             # project folder on Sherlock
OUT <- file.path(W, "output")
dir.create(file.path(OUT, "loop"), showWarnings = FALSE)
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))

# log_(): print a time-stamped message
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")


# ---- Settings -------------------------------------------------------------------------------------------------
## MIN_DROP: author 2026-10-14 first 200 (r3 -> r7), then "until < 100" (continued from r7 with START_ROUND=7 MIN_DROP=100)
MIN_DROP <- as.integer(Sys.getenv("MIN_DROP", "200"))
MAX_ROUNDS <- 10                                                     # at most 10 rounds per run
FROM <- as.integer(Sys.getenv("START_ROUND", "3"))                   # start round (3 = the r3 object)

## ROLLING=1 (author 2026-10-14): after each round every kept cell is re-annotated with its new cluster's type, so the next
## round compares against the previous round's annotation (protected HRS / T_gd cells are re-typed too; they stay protected)
## [note added in the 2026-10-23 rewrite: the code in step 2f does NOT re-type protected cells (roll <- !keep), as the
##  protected-types note below says; the "re-typed too" remark above is older than that code]
ROLL <- Sys.getenv("ROLLING", "0") == "1"

# file-name suffix so that each run keeps its own tables: "" (default run), "_from_r7", "_rolling", "_from_r7_rolling"
from_suffix <- if (FROM == 3L) "" else paste0("_from_r", FROM)
rolling_suffix <- if (ROLL) "_rolling" else ""
SUF <- paste0(from_suffix, rolling_suffix)

CMP <- "cell_type"                                    # comparison level ("lineage" = the superseded first run)

## r3 cluster -> cell type (majority published label, cpA3 composition table 2026-10-14) -> lineage ----------------
TYPE_BY_R3 <- c(`0` = "B_preGC", `1` = "B_preGC", `5` = "B_Memory", `14` = "B_GC", `15` = "B_GC", `12` = "Plasma", `18` = "Plasma",
                `2` = "T_CD4", `6` = "T_CD4", `3` = "T_Treg", `4` = "T_CD8", `8` = "NK", `9` = "T_MKI67",
                `7` = "Macrophage", `11` = "DC_Conventional", `16` = "DC_Conventional", `10` = "DC_Plasmacytoid",
                `13` = "CAF", `17` = "HRS")
LINEAGE <- c(T_CD4 = "T_NK", T_Treg = "T_NK", T_CD8 = "T_NK", NK = "T_NK", T_gd = "T_NK", ILC = "T_NK", T_MKI67 = "T_NK",
             B_preGC = "B", B_Memory = "B", B_GC = "B", Plasma = "B",
             Macrophage = "Myeloid", DC_Conventional = "Myeloid", DC_Plasmacytoid = "Myeloid", CAF = "Stromal", HRS = "HRS")

## protected rare types (author 2026-10-14): HRS, T_gd, ILC (T_gd never forms its own cluster; ILC has none at r3). Protected
## cells (by their own published label) are never removed, start the loop with their own label (not their r3 cluster's type)
## and are never re-typed by the rolling step, as LP / Treg in the NLPHL loop.
ALWAYS_KEEP_TYPES <- "HRS"                            # every cell of the HRS cluster type
ALWAYS_KEEP_LABELS <- c("HRS", "T_gd", "ILC")         # every cell published with one of these labels


# ---- recluster(): the rounds 1-3 recipe on the current cells, with a given seed ------------------------------------
#   per-dataset 3,000 HVGs ranked by how many datasets pick them -> ScaleData -> PCA 50 (seed) -> Harmony (dataset +
#   chemistry; set.seed(seed)) -> neighbours on 30 dims -> Louvain res 0.6 (seed), or res 0.3 / 0.6 / 1.0 when final = TRUE.
#   The chosen HVGs are attached to the returned object as attribute "hvg".
recluster <- function(obj, seed, final = FALSE) {
  obj@reductions <- list()
  obj@graphs <- list()
  obj@neighbors <- list()

  # 3,000 HVGs within each dataset (datasets in alphabetical order)
  cells_by_dataset <- split(colnames(obj), obj$dataset)
  hvg_by_dataset <- list()
  for (dataset in names(cells_by_dataset)) {
    dataset_counts <- LayerData(obj, assay = "RNA", layer = "counts")[, cells_by_dataset[[dataset]]]
    dataset_obj <- FindVariableFeatures(CreateSeuratObject(dataset_counts), nfeatures = 3000, verbose = FALSE)
    hvg_by_dataset[[dataset]] <- VariableFeatures(dataset_obj)
  }

  # for every gene picked anywhere: in how many datasets it was picked, and its median rank in those datasets
  n_datasets_per_gene <- table(unlist(hvg_by_dataset))
  candidate_genes <- names(n_datasets_per_gene)
  median_rank <- numeric(length(candidate_genes))
  for (i in seq_along(candidate_genes)) {
    ranks <- c()
    for (dataset in names(hvg_by_dataset)) {
      ranks <- c(ranks, match(candidate_genes[i], hvg_by_dataset[[dataset]]))   # NA when not picked in this dataset
    }
    median_rank[i] <- median(ranks, na.rm = TRUE)
  }

  # top 3,000: most datasets first, ties broken by the better (lower) median rank
  hvg <- candidate_genes[order(-n_datasets_per_gene, median_rank)][1:3000]
  DefaultAssay(obj) <- "RNA"
  VariableFeatures(obj) <- hvg

  obj <- ScaleData(obj, features = hvg, verbose = FALSE)
  obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = seed, verbose = FALSE)
  set.seed(seed)                                                       # seed for Harmony
  obj <- RunHarmony(obj, group.by.vars = c("dataset", "chemistry"), reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
  resolutions <- if (final) c(0.3, 0.6, 1.0) else 0.6
  for (resolution in resolutions) {
    obj <- FindClusters(obj, resolution = resolution, random.seed = seed, verbose = FALSE)
  }
  attr(obj, "hvg") <- hvg
  obj
}


# ---- Step 1. Starting object and starting cell types ----------------------------------------------------------------
if (FROM == 3L) {
  obj <- readRDS(file.path(OUT, "cHL_merged_harmony_v2_r3.rds"))
  obj$r3_cluster <- as.character(obj$RNA_snn_res.0.6)
  stopifnot(setequal(names(TYPE_BY_R3), unique(obj$r3_cluster)))     # every r3 cluster has a cell type
  # every cell starts with its r3 cluster's cell type
  obj$cell_type <- unname(TYPE_BY_R3[obj$r3_cluster])
  obj$lineage <- unname(LINEAGE[obj$cell_type])
  stopifnot(!anyNA(obj$lineage))
  # protected cells: never removed (loop_keep), and they start with their own published label
  obj$loop_keep <- obj$cell_type %in% ALWAYS_KEEP_TYPES | obj$label %in% ALWAYS_KEEP_LABELS
  own <- obj$label %in% ALWAYS_KEEP_LABELS
  obj$cell_type[own] <- obj$label[own]                                # protected labels kept as they are
  obj$lineage <- unname(LINEAGE[obj$cell_type])
} else {                                                     # continue: the r3-derived cell types stay fixed, as in the ex vivo loop
  obj <- readRDS(file.path(OUT, sprintf("cHL_merged_harmony_v2_r%d.rds", FROM)))
  stopifnot(!anyNA(obj$cell_type), !anyNA(obj$loop_keep))
}
log_("start r%d: %d cells; always kept: %d; stop when a round drops < %d", FROM, ncol(obj), sum(obj$loop_keep), MIN_DROP)


# ---- Step 2. The loop: recluster, find the unstable cells, drop them --------------------------------------------------
rounds <- list()          # one summary row per round
removed <- list()         # the removed cells of each round
n <- FROM                 # round number
repeat {
  n <- n + 1L
  obj <- recluster(obj, seed = n)                              # a new draw (r3 itself used seeds 42 / 0)

  # 2a. per-cell table for this round: new cluster (cl), the annotation compared (cmp), protected (keep)
  md <- data.table(cell = colnames(obj), cl = as.character(obj$RNA_snn_res.0.6), lineage = obj$lineage,
                   cmp = obj@meta.data[[CMP]], keep = obj$loop_keep,
                   dataset = obj$dataset, sample = obj$sample, label = obj$label, cell_type = obj$cell_type)

  # 2b. type of each new cluster = the most common annotation among its cells (and its %)
  type_counts <- md[, .N, by = .(cl, cmp)]
  type_counts <- type_counts[order(cl, -N)]
  ct <- type_counts[, .(cluster_type = cmp[1], type_pct = round(100 * N[1] / sum(N), 1)), by = cl]
  md <- merge(md, ct, by = "cl")

  # 2c. unstable = annotation differs from its cluster's type, and the cell is not protected
  md[, unstable := cmp != cluster_type & !keep]
  n_drop <- sum(md$unstable)
  n_protected_mismatch <- md[cmp != cluster_type & keep, .N]
  log_("round r%d: %d cells, %d clusters, %d unstable (+%d protected cells would have been unstable)",
       n, nrow(md), nrow(ct), n_drop, n_protected_mismatch)
  unstable_moves <- md[unstable == TRUE, .N, by = .(cell_type, cluster_type)]
  print(unstable_moves[order(-N)])

  # 2d. round summary: cells in, clusters, dropped, protected mismatches, % dropped per dataset
  pct <- md[, .(pct = round(100 * mean(unstable), 2)), by = dataset]
  pct <- pct[order(dataset)]
  round_summary <- list(round = paste0("r", n), cells_in = nrow(md), clusters = nrow(ct), dropped = n_drop,
                        protected_mismatch = n_protected_mismatch)
  for (i in seq_len(nrow(pct))) {
    round_summary[[paste0("pct_drop_", pct$dataset[i])]] <- pct$pct[i]
  }
  rounds[[n]] <- as.data.table(round_summary)
  print(rbindlist(rounds, fill = TRUE))

  # 2e. record and drop the unstable cells
  removed[[n]] <- md[unstable == TRUE, .(cell, round = paste0("r", n), dataset, sample, label, cell_type, lineage,
                                         new_cluster = cl, cluster_type)]
  obj <- subset(obj, cells = md[unstable == FALSE, cell])

  # 2f. rolling (ROLLING=1 only): kept, non-protected cells take their new cluster's type
  if (ROLL) {
    kept_rows <- md[match(colnames(obj), cell)]
    roll <- !kept_rows$keep                                           # protected cells keep their own label
    obj$cell_type[roll] <- kept_rows$cluster_type[roll]
    obj$lineage <- unname(LINEAGE[obj$cell_type])
    log_("rolling: %d non-protected cells re-typed (should be 0: kept cells already match their cluster)",
         kept_rows[roll & cmp != cluster_type, .N])
  }

  # 2g. save the kept cell ids of this round; stop when few cells were dropped or after MAX_ROUNDS rounds
  saveRDS(data.table(cell = colnames(obj)), file.path(OUT, "loop", sprintf("cell_ids_r%d%s.rds", n, SUF)))
  if (n_drop < MIN_DROP || n - FROM >= MAX_ROUNDS) break
}
if (n_drop >= MIN_DROP) cat("WARNING: stopped at MAX_ROUNDS before a round removed fewer than", MIN_DROP, "cells\n")
fwrite(rbindlist(rounds, fill = TRUE), file.path(OUT, paste0("cpA_loop_rounds", SUF, ".tsv")), sep = "\t")
fwrite(rbindlist(removed), file.path(OUT, paste0("cpA_loop_removed_cells", SUF, ".tsv.gz")), sep = "\t")


# ---- Step 3. Final round rN: full recluster on the kept cells + UMAP + exports (same layout as cpA3_) ---------------------
TAGN <- paste0("r", n, if (ROLL) "_rolling" else "")         # object name tag, e.g. "r8" or "r8_rolling"
PFX <- paste0("cpA", n, if (ROLL) "rolling" else "")         # table prefix, e.g. "cpA8" or "cpA8rolling"

# 3a. recluster with seed = round number + 100, at res 0.3 / 0.6 / 1.0; save the HVGs
obj <- recluster(obj, seed = n + 100L, final = TRUE)
hvg <- attr(obj, "hvg")
fwrite(data.table(ensg = hvg, symbol = v47$symbol[match(hvg, v47$ensg)]), file.path(OUT, paste0(PFX, "_hvg_3000.tsv")), sep = "\t")

# 3b. UMAPs (seed 42): unintegrated (PCA) and Harmony
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 42, reduction.name = "umap_unintegrated", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)

# 3c. mixing: neighbour diversity (inverse Simpson over k = 30 harmony neighbours)
#     For every cell: 1 / sum(p^2), where p = the share of its 30 nearest neighbours from each dataset (or chemistry).
neighbours <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, k.param = 30, return.neighbor = TRUE, verbose = FALSE)@neighbors[[1]]
neighbour_index <- Indices(neighbours)                  # cells x 30 matrix of neighbour positions
inverse_simpson <- function(groups) {
  scores <- numeric(nrow(neighbour_index))
  for (i in seq_len(nrow(neighbour_index))) {
    neighbour_cells <- neighbour_index[i, ]
    proportions <- table(groups[neighbour_cells]) / length(neighbour_cells)
    scores[i] <- 1 / sum(proportions^2)
  }
  scores
}
obj$ilisi_dataset <- inverse_simpson(obj$dataset)
obj$ilisi_chemistry <- inverse_simpson(obj$chemistry)
log_("final round %s: %d cells", TAGN, ncol(obj))

# 3d. per-cell table: metadata, earlier clusters, loop cell type / lineage, QC, all cluster columns, mixing, UMAP coordinates
umap_harmony <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
cluster_cols <- grep("^RNA_snn_res", names(obj@meta.data), value = TRUE)
cells <- data.table(obj@meta.data[, c("dataset", "sample", "chemistry", "label", "label_source", "s1_cell_type", "s1_cell_state",
                                      "s1_state_matches_type", "ref_exclude", "r1_cluster_res0.6", "r2_cluster_res0.6", "major_type",
                                      "r2_subcluster", "r3_cluster", "cell_type", "lineage", cluster_cols, "ilisi_dataset",
                                      "ilisi_chemistry", "nCount_RNA", "nFeature_RNA", "nCount_orig", "nFeature_orig", "percent.mt",
                                      "contaminated", "contam_lineages", "protected")],
                    keep.rownames = "cell")
cells[, UMAP_1 := umap_harmony[, 1]]
cells[, UMAP_2 := umap_harmony[, 2]]
cells[, UMAPu_1 := umap_unintegrated[, 1]]
cells[, UMAPu_2 := umap_unintegrated[, 2]]
fwrite(cells, file.path(OUT, paste0(PFX, "_cells.tsv.gz")), sep = "\t")

# 3e. composition check: for each final cluster (res 0.6), cells and % per loop cell type
chk <- cells[, .N, by = .(cluster = RNA_snn_res.0.6, cell_type)]
chk[, pct := round(100 * N / sum(N), 1), by = cluster]
chk <- chk[order(cluster, -N)]
fwrite(chk, file.path(OUT, paste0(PFX, "_cluster_by_cell_type.tsv")), sep = "\t")
# print each cluster's main cell type, least pure clusters first
main_type <- chk[, .SD[1], by = cluster]
print(main_type[order(pct)])

# 3f. dot-plot tables: for each marker, mean log-expression and % of cells expressing it, per group
markers <- c("PTPRC", "CD3E", "CD4", "CD8A", "FOXP3", "IL2RA", "MKI67", "TOP2A", "NKG7", "GNLY", "KLRD1", "TRDC", "TRGC1", "KIT",
             "IL7R", "MS4A1", "CD19", "TCL1A", "IGHD", "CD27", "BCL6", "AICDA", "RGS13", "MZB1", "JCHAIN", "SDC1", "CD68", "CD163",
             "C1QA", "LYZ", "CD14", "VCAN", "FCGR3A", "CLEC9A", "CD1C", "LILRA4", "IL3RA", "COL1A1", "COL1A2", "PDGFRA", "PECAM1",
             "VWF", "TNFRSF8", "IRF4", "CD274", "FOS", "JUN", "HSPA1A", "EGR1")
marker_genes <- v47[symbol %in% markers & ensg %in% rownames(obj)]          # in gene-table order
marker_data <- LayerData(obj, assay = "RNA", layer = "data")[marker_genes$ensg, ]   # log-normalised values

# dotplot_table(): one block of rows (one per marker) per group; column "group" = the group name
dotplot_table <- function(groups) {
  cell_positions <- split(seq_len(ncol(obj)), groups)
  blocks <- list()
  for (group_name in names(cell_positions)) {
    group_data <- marker_data[, cell_positions[[group_name]], drop = FALSE]
    blocks[[group_name]] <- data.table(gene = marker_genes$symbol,
                                       mean_logexp = Matrix::rowMeans(group_data),
                                       pct_expr = 100 * Matrix::rowMeans(group_data > 0))
  }
  rbindlist(blocks, idcol = "group")
}
fwrite(dotplot_table(obj$RNA_snn_res.0.6), file.path(OUT, paste0(PFX, "_dotplot_by_cluster_res0.6.tsv")), sep = "\t")
# by published label; S1 cells without a label are grouped by their provisional type ("S1:<type>")
label_groups <- ifelse(is.na(obj$label), paste0("S1:", obj$s1_cell_type), obj$label)
fwrite(dotplot_table(label_groups), file.path(OUT, paste0(PFX, "_dotplot_by_label.tsv")), sep = "\t")

# 3g. save the final object
saveRDS(obj, file.path(OUT, paste0("cHL_merged_harmony_v2_", TAGN, ".rds")))
log_("saved %s object and %s_ tables", TAGN, PFX)
