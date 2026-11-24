# =============================================================================
# v2_C4a_s1_cluster_check.R
#
# PURPOSE : Quick check of four round-10 clusters (res 0.6 clusters 10, 14, 9 and 8): for the main groups of cells in each
#           cluster (S1 cells by their ex vivo cell type, other cells by their published label), print the % of cells
#           that express each of a set of marker genes (raw count > 0).
# INPUT   : output/cHL_merged_harmony_v2_r10.rds, ref/gene_table_v47.tsv (on OAK)
# OUTPUT  : printed tables only (no files written)
# RUN     : on Sherlock, Rscript v2_C4a_s1_cluster_check.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})
W <- "${PROJECT_ROOT}/cHL_v2"      # project root on OAK
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))
obj <- readRDS(file.path(W, "output/cHL_merged_harmony_v2_r10.rds"))

# Step 1. Raw counts of the marker genes (object rows are Ensembl IDs; genes absent from the object are skipped)
marker_genes <- c("MKI67", "TOP2A", "CD3E", "CD4", "CD8A", "FOXP3", "KIT", "IL7R", "NKG7", "GNLY", "KLRF1", "TRDC", "IL4I1",
                  "RORC", "CLEC4C", "LILRA4", "CD1C", "TNFRSF8", "CCL17")
marker_ids <- v47$ensg[match(marker_genes, v47$symbol)]
in_object <- marker_ids %in% rownames(obj)
marker_counts <- LayerData(obj, assay = "RNA", layer = "counts")[marker_ids[in_object], ]
rownames(marker_counts) <- marker_genes[in_object]

# Step 2. One row per cell: dataset, S1 ex vivo type, published label, round-10 cluster (res 0.6)
cell_table <- data.table(obj@meta.data[, c("dataset", "s1_cell_type", "label")], keep.rownames = "cell")
cell_table[, cl := as.character(obj$RNA_snn_res.0.6)]

# Step 3. For each cluster: the (up to) 7 largest groups with >= 20 cells, and % of their cells expressing each marker
for (k in c("10", "14", "9", "8")) {
  # group = "S1:<ex vivo type>" for S1 cells, "nonS1:<published label>" for the others
  cluster_cells <- cell_table[cl == k, .(cell, g = fifelse(dataset == "s1", paste0("S1:", s1_cell_type), paste0("nonS1:", label)))]
  group_sizes <- cluster_cells[, .N, by = g]
  group_sizes <- group_sizes[N >= 20]
  group_sizes <- group_sizes[order(-N)]
  top <- group_sizes[1:min(.N, 7)]                    # largest 7 groups

  pct_table <- matrix(NA_real_, nrow = nrow(marker_counts), ncol = nrow(top),
                      dimnames = list(rownames(marker_counts), NULL))
  for (col in seq_len(nrow(top))) {
    group_cells <- cluster_cells[g == top$g[col], cell]
    pct_table[, col] <- round(100 * Matrix::rowMeans(marker_counts[, group_cells, drop = FALSE] > 0))
  }
  colnames(pct_table) <- paste0(top$g, "(", top$N, ")")    # column name = group (number of cells)
  cat("\n== r10 cluster", k, ": % cells expressing\n")
  print(pct_table)
}
