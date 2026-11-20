# =============================================================================
# v2_C3_recluster_r4.R — round r4 (Sherlock): remove the S1 cells with only a provisional label and re-integrate.
#
# Author decision 2026-10-14: remove the S1 cells that carry only a provisional label
#   (s1_cell_type Endothelial, IEG_stress, Innate_T, Monocyte; label = NA), then redo
#   HVG -> PCA -> Harmony (dataset + chemistry) -> clusters -> UMAP (round-1 recipe).
# Recipe = rounds 1-3: 3,000 HVGs chosen per dataset and ranked by how many datasets pick them -> ScaleData -> PCA 50 ->
#   Harmony (dataset + chemistry) -> neighbours on 30 dims -> Louvain res 0.3 / 0.6 / 1.0 -> UMAP.
#   Seeds: 42 (PCA, UMAP, Harmony), 0 (clustering).
# Input : output/cHL_merged_harmony_v2_r3.rds
# Output: output/cHL_merged_harmony_v2_r4.rds, output/cpA4_*.tsv(.gz), output/cpA4_removed_cells.tsv.gz
# Run   : Rscript v2_C3_recluster_r4.R
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

# S1 provisional types to remove (author 2026-10-14)
PROV <- c("Endothelial", "IEG_stress", "Innate_T", "Monocyte")


# ---- Step 1. Load r3 and list the cells to remove --------------------------------------------------------------
obj <- readRDS(file.path(OUT, "cHL_merged_harmony_v2_r3.rds"))
cell_meta <- obj@meta.data
is_provisional <- is.na(cell_meta$label) & cell_meta$dataset == "s1" & cell_meta$s1_cell_type %in% PROV
rm_cells <- rownames(cell_meta)[is_provisional]
stopifnot(all(cell_meta$s1_cell_type[is.na(cell_meta$label)] %in% PROV))          # no other unlabelled cells

# save the removed cells: cell, reason, S1 type, sample, r3 cluster
rm_tab <- data.table(cell = rm_cells, reason = "S1 provisional label (author 2026-10-14)",
                     s1_cell_type = cell_meta[rm_cells, "s1_cell_type"],
                     sample = cell_meta[rm_cells, "sample"],
                     r3_cluster = as.character(cell_meta[rm_cells, "RNA_snn_res.0.6"]))
rm(cell_meta)
fwrite(rm_tab, file.path(OUT, "cpA4_removed_cells.tsv.gz"), sep = "\t")
print(rm_tab[, .N, by = s1_cell_type])


# ---- Step 2. Drop them; keep the r3 cluster as r3_cluster_res0.6; clear old clusters and reductions -----------------
obj <- subset(obj, cells = setdiff(colnames(obj), rm_cells))
obj$r3_cluster_res0.6 <- obj$RNA_snn_res.0.6
old_cluster_cols <- grep("^RNA_snn_res|^seurat_clusters$", names(obj@meta.data), value = TRUE)
for (col in old_cluster_cols) {
  obj[[col]] <- NULL
}
obj@reductions <- list()
obj@graphs <- list()
obj@neighbors <- list()
log_("removed %d S1 provisional cells; %d cells remain", length(rm_cells), ncol(obj))


# ---- Step 3. Choose 3,000 HVGs (same integration recipe as rounds 1-3) -------------------------------------------------
# 3a. 3,000 HVGs within each dataset (datasets in alphabetical order)
cells_by_dataset <- split(colnames(obj), obj$dataset)
hvg_by_dataset <- list()
for (dataset in names(cells_by_dataset)) {
  dataset_counts <- LayerData(obj, assay = "RNA", layer = "counts")[, cells_by_dataset[[dataset]]]
  dataset_obj <- FindVariableFeatures(CreateSeuratObject(dataset_counts), nfeatures = 3000, verbose = FALSE)
  hvg_by_dataset[[dataset]] <- VariableFeatures(dataset_obj)
}

# 3b. For every gene picked anywhere: in how many datasets it was picked, and its median rank in those datasets
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

# 3c. Top 3,000: most datasets first, ties broken by the better (lower) median rank
hvg <- candidate_genes[order(-n_datasets_per_gene, median_rank)][1:3000]
DefaultAssay(obj) <- "RNA"
VariableFeatures(obj) <- hvg
fwrite(data.table(ensg = hvg, symbol = v47$symbol[match(hvg, v47$ensg)], n_datasets = as.integer(n_datasets_per_gene[hvg])),
       file.path(OUT, "cpA4_hvg_3000.tsv"), sep = "\t")


# ---- Step 4. Scale, PCA, Harmony, clusters, UMAP -------------------------------------------------------------------
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


# ---- Step 5. Mixing: neighbour diversity (inverse Simpson over k = 30 harmony neighbours) ---------------------------
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
log_("round-4 integration done")


# ---- Step 6. Exports (same layout as rounds 1-3, prefix cpA4_) + r4 cluster x r2 major type check -----------------------
# 6a. per-cell table: metadata, earlier clusters, QC, all cluster columns, mixing scores, UMAP coordinates
umap_harmony <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
cluster_cols <- grep("^RNA_snn_res", names(obj@meta.data), value = TRUE)
cells <- data.table(obj@meta.data[, c("dataset", "sample", "chemistry", "label", "label_source", "s1_cell_type", "s1_cell_state",
                                      "s1_state_matches_type", "ref_exclude", "r1_cluster_res0.6", "r2_cluster_res0.6",
                                      "r3_cluster_res0.6", "major_type", "r2_subcluster", cluster_cols, "ilisi_dataset",
                                      "ilisi_chemistry", "nCount_RNA", "nFeature_RNA", "nCount_orig", "nFeature_orig",
                                      "percent.mt", "contaminated", "contam_lineages", "protected")],
                    keep.rownames = "cell")
cells[, UMAP_1 := umap_harmony[, 1]]
cells[, UMAP_2 := umap_harmony[, 2]]
cells[, UMAPu_1 := umap_unintegrated[, 1]]
cells[, UMAPu_2 := umap_unintegrated[, 2]]
fwrite(cells, file.path(OUT, "cpA4_cells.tsv.gz"), sep = "\t")

# 6b. composition check: for each new r4 cluster (res 0.6), cells and % per r2 major type
major_by_cluster <- cells[, .N, by = .(r4_cluster = RNA_snn_res.0.6, major_type)]
major_by_cluster[, pct := round(100 * N / sum(N), 1), by = r4_cluster]
major_by_cluster <- major_by_cluster[order(r4_cluster, -N)]
fwrite(major_by_cluster, file.path(OUT, "cpA4_r4cluster_by_r2_major_type.tsv"), sep = "\t")
# print each cluster's main major type, least pure clusters first
main_major <- major_by_cluster[, .SD[1], by = r4_cluster]
print(main_major[order(pct)])

# 6c. dot-plot tables: for each marker, mean log-expression and % of cells expressing it, per group
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
fwrite(dotplot_table(obj$RNA_snn_res.0.6), file.path(OUT, "cpA4_dotplot_by_cluster_res0.6.tsv"), sep = "\t")
# by published label; S1 cells without a label are grouped by their provisional type ("S1:<type>")
label_groups <- ifelse(is.na(obj$label), paste0("S1:", obj$s1_cell_type), obj$label)
fwrite(dotplot_table(label_groups), file.path(OUT, "cpA4_dotplot_by_label.tsv"), sep = "\t")

# 6d. save the round-4 object
saveRDS(obj, file.path(OUT, "cHL_merged_harmony_v2_r4.rds"))
log_("saved round-4 object and tables")
