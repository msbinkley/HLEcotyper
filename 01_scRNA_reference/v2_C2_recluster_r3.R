# =============================================================================
# v2_C2_recluster_r3.R — round r3 (Sherlock): drop the trash sub-UMAP subclusters and re-integrate.
#
# Author decisions 2026-10-14 (figures_v2/subUMAP/filter_decisions_v2.md):
#   drop the r2 sub-UMAP subclusters listed in DROP below (A: trash, B: heat-shock, C: HRS + monocyte subcluster), then redo
#   HVG -> PCA -> Harmony (dataset + chemistry) -> clusters -> UMAP on the remaining cells (round-1 recipe).
# Recipe = rounds 1-2: 3,000 HVGs chosen per dataset and ranked by how many datasets pick them -> ScaleData -> PCA 50 ->
#   Harmony (dataset + chemistry) -> neighbours on 30 dims -> Louvain res 0.3 / 0.6 / 1.0 -> UMAP.
#   Seeds: 42 (PCA, UMAP, Harmony), 0 (clustering).
# Input : output/cHL_merged_v2_r2_qc.rds (r2 + QC/lineage scores/major type), output/subumap/r2_sub_<TYPE>_cell_ids.rds
# Output: output/cHL_merged_harmony_v2_r3.rds, output/cpA3_*.tsv(.gz), output/cpA3_removed_cells.tsv.gz
# Run   : Rscript v2_C2_recluster_r3.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
})
set.seed(42)                                                          # global seed
W <- "${PROJECT_ROOT}/cHL_v2"             # project folder on Sherlock
OUT <- file.path(W, "output")
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))

# log_(): print a time-stamped message
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

# sub-UMAP subclusters to drop: major type, subcluster, reason (author 2026-10-14)
DROP <- rbindlist(list(
  list("Stromal",        "2",  "A: CAF-macrophage doublets (macrophage genes, myeloid score high)"),
  list("Stromal",        "5",  "A: mixed CAF/T labels, T genes"),
  list("T_MKI67",        "10", "A: low quality (MT- genes only, 13 significant markers)"),
  list("T_MKI67",        "11", "A: B_GC cells inside T_MKI67"),
  list("T_CD8_NK_gd",    "13", "A: low quality (metallothionein/MT-/JUN, 42 significant markers)"),
  list("T_CD4_Treg",     "4",  "B: heat-shock (dissociation stress)"),
  list("B_preGC_Memory", "9",  "B: heat-shock (dissociation stress)"),
  list("T_CD8_NK_gd",    "6",  "B: heat-shock (dissociation stress)"),
  list("Plasma_BGC_HRS", "9",  "C: HRS with monocyte program (HRS-monocyte doublets/rosettes)")))
setnames(DROP, c("major_type", "subcluster", "reason"))


# ---- Step 1. Load r2 and every cell's sub-UMAP subcluster -------------------------------------------------------
obj <- readRDS(file.path(OUT, "cHL_merged_v2_r2_qc.rds"))

# one cell-id table per major type (written by v2_C1_subumap.R), stacked, with the major type added
subcluster_tables <- list()
for (major in unique(obj$major_type)) {
  type_ids <- readRDS(file.path(OUT, "subumap", paste0("r2_sub_", major, "_cell_ids.rds")))
  type_ids[, major_type := major]
  subcluster_tables[[length(subcluster_tables) + 1]] <- type_ids
}
sub <- rbindlist(subcluster_tables)
stopifnot(setequal(sub$cell, colnames(obj)), !anyDuplicated(sub$cell))   # every cell in exactly one sub-UMAP
obj$r2_subcluster <- paste0(sub$major_type, ":", sub$subcluster)[match(colnames(obj), sub$cell)]   # e.g. "Stromal:2"


# ---- Step 2. The cells to remove (cells of the DROP subclusters) -----------------------------------------------------
rm_tab <- merge(sub, DROP, by = c("major_type", "subcluster"))
rm_cells <- rm_tab$cell
cell_meta <- obj@meta.data
rm_tab[, dataset := cell_meta[cell, "dataset"]]
rm_tab[, sample := cell_meta[cell, "sample"]]
rm_tab[, label := cell_meta[cell, "label"]]
rm(cell_meta)
fwrite(rm_tab, file.path(OUT, "cpA3_removed_cells.tsv.gz"), sep = "\t")
print(rm_tab[, .N, by = .(major_type, subcluster, reason)])
print(rm_tab[label %in% c("T_gd", "HRS"), .N, by = .(major_type, subcluster, label)])   # protected cells among them


# ---- Step 3. Drop them; keep the r2 cluster as r2_cluster_res0.6; clear old clusters and reductions -----------------
obj <- subset(obj, cells = setdiff(colnames(obj), rm_cells))
obj$r2_cluster_res0.6 <- obj$RNA_snn_res.0.6
old_cluster_cols <- grep("^RNA_snn_res|^seurat_clusters$", names(obj@meta.data), value = TRUE)
for (col in old_cluster_cols) {
  obj[[col]] <- NULL
}
obj@reductions <- list()
obj@graphs <- list()
obj@neighbors <- list()
log_("removed %d cells (%d subclusters); %d cells remain", length(rm_cells), nrow(DROP), ncol(obj))


# ---- Step 4. Choose 3,000 HVGs (same integration recipe as rounds 1-2) -------------------------------------------------
# 4a. 3,000 HVGs within each dataset (datasets in alphabetical order)
cells_by_dataset <- split(colnames(obj), obj$dataset)
hvg_by_dataset <- list()
for (dataset in names(cells_by_dataset)) {
  dataset_counts <- LayerData(obj, assay = "RNA", layer = "counts")[, cells_by_dataset[[dataset]]]
  dataset_obj <- FindVariableFeatures(CreateSeuratObject(dataset_counts), nfeatures = 3000, verbose = FALSE)
  hvg_by_dataset[[dataset]] <- VariableFeatures(dataset_obj)
}

# 4b. For every gene picked anywhere: in how many datasets it was picked, and its median rank in those datasets
n_datasets_per_gene <- table(unlist(hvg_by_dataset))
candidate_genes <- names(n_datasets_per_gene)
median_rank <- numeric(length(candidate_genes))
for (i in seq_along(candidate_genes)) {
  ranks <- c()
  for (dataset in names(hvg_by_dataset)) {
    ranks <- c(ranks, match(candidate_genes[i], hvg_by_dataset[[dataset]]))     # NA when not picked in this dataset
  }
  median_rank[i] <- median(ranks, na.rm = TRUE)
}

# 4c. Top 3,000: most datasets first, ties broken by the better (lower) median rank
hvg <- candidate_genes[order(-n_datasets_per_gene, median_rank)][1:3000]
DefaultAssay(obj) <- "RNA"
VariableFeatures(obj) <- hvg
fwrite(data.table(ensg = hvg, symbol = v47$symbol[match(hvg, v47$ensg)], n_datasets = as.integer(n_datasets_per_gene[hvg])),
       file.path(OUT, "cpA3_hvg_3000.tsv"), sep = "\t")


# ---- Step 5. Scale, PCA, Harmony, clusters, UMAP -------------------------------------------------------------------
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = 42, verbose = FALSE)                          # 50 PCs, seed 42
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 42, reduction.name = "umap_unintegrated", verbose = FALSE)
set.seed(42)                                                                                           # seed for Harmony
obj <- RunHarmony(obj, group.by.vars = c("dataset", "chemistry"), reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)                         # 30 Harmony dims
for (resolution in c(0.3, 0.6, 1.0)) {                                                                 # Louvain (default), seed 0
  obj <- FindClusters(obj, resolution = resolution, random.seed = 0, verbose = FALSE)
}
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)


# ---- Step 6. Mixing: neighbour diversity (inverse Simpson over k = 30 harmony neighbours) ---------------------------
# For every cell: 1 / sum(p^2), where p = the share of its 30 nearest neighbours (Harmony space) from each dataset
# (or chemistry). 1 = all neighbours from one group; higher = better mixed.
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
log_("round-3 integration done")


# ---- Step 7. Exports (same layout as rounds 1-2, prefix cpA3_) + r3 cluster x r2 major type check -----------------------
# 7a. per-cell table: metadata, earlier clusters, QC, all cluster columns, mixing scores, UMAP coordinates
umap_harmony <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
cluster_cols <- grep("^RNA_snn_res", names(obj@meta.data), value = TRUE)
cells <- data.table(obj@meta.data[, c("dataset", "sample", "chemistry", "label", "label_source", "s1_cell_type", "s1_cell_state",
                                      "s1_state_matches_type", "ref_exclude", "r1_cluster_res0.6", "r2_cluster_res0.6", "major_type",
                                      "r2_subcluster", cluster_cols, "ilisi_dataset", "ilisi_chemistry", "nCount_RNA", "nFeature_RNA",
                                      "nCount_orig", "nFeature_orig", "percent.mt", "contaminated", "contam_lineages", "protected")],
                    keep.rownames = "cell")
cells[, UMAP_1 := umap_harmony[, 1]]
cells[, UMAP_2 := umap_harmony[, 2]]
cells[, UMAPu_1 := umap_unintegrated[, 1]]
cells[, UMAPu_2 := umap_unintegrated[, 2]]
fwrite(cells, file.path(OUT, "cpA3_cells.tsv.gz"), sep = "\t")

# 7b. composition check: for each new r3 cluster (res 0.6), cells and % per r2 major type
major_by_cluster <- cells[, .N, by = .(r3_cluster = RNA_snn_res.0.6, major_type)]
major_by_cluster[, pct := round(100 * N / sum(N), 1), by = r3_cluster]
major_by_cluster <- major_by_cluster[order(r3_cluster, -N)]
fwrite(major_by_cluster, file.path(OUT, "cpA3_r3cluster_by_r2_major_type.tsv"), sep = "\t")
# print each cluster's main major type, least pure clusters first
main_major <- major_by_cluster[, .SD[1], by = r3_cluster]
print(main_major[order(pct)])

# 7c. dot-plot tables: for each marker, mean log-expression and % of cells expressing it, per group
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
fwrite(dotplot_table(obj$RNA_snn_res.0.6), file.path(OUT, "cpA3_dotplot_by_cluster_res0.6.tsv"), sep = "\t")
# by published label; S1 cells without a label are grouped by their provisional type ("S1:<type>")
label_groups <- ifelse(is.na(obj$label), paste0("S1:", obj$s1_cell_type), obj$label)
fwrite(dotplot_table(label_groups), file.path(OUT, "cpA3_dotplot_by_label.tsv"), sep = "\t")

# 7d. save the round-3 object
saveRDS(obj, file.path(OUT, "cHL_merged_harmony_v2_r3.rds"))
log_("saved round-3 object and tables")
