# =============================================================================
# v2_G5_exvivo_B_subcluster.R — sub-cluster the ex vivo "B" cells into B_preGC / B_postGC before the EcoTyper
# recovery of the final round-14 model (Fig 8e); author 2026-10-24: "Sub-cluster B first".
#
# Why: the ex vivo labels (f6 object, column cell_type) have one "B" type, but the final model has B_preGC and
#   B_postGC (= B_GC + B_Memory) states.
#
# Recipe (author's sub-UMAP rules; one dataset / one chemistry, so LogNormalize rather than SCTransform):
#   B cells of the f6 object -> LogNormalize -> 2,000 HVG (vst) -> scale -> PCA 30 -> Harmony on sample
#   -> Louvain res 0.4 -> UMAP. Seed 42 everywhere.
#   Markers by the author rule: FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) on the
#   log-normalised RNA, ranked by p_adj then avg_log2FC, p_adj < 0.05, top 50 per cluster (+ how many of 50 pass).
#   preGC / postGC scores: Seurat module scores of marker genes taken from the ROUND-14 CIBERSORTx signature itself
#   (sigmtx_merged.txt): the 30 signature genes with the largest log2 ratio B_preGC / mean(B_GC, B_Memory) = "preGC"
#   genes, and the 30 with the largest inverse ratio = "postGC" genes (only genes present in the ex vivo object).
#   PROPOSED label per cluster = the higher mean score. Nothing is used until the author has seen the UMAP, the
#   markers and the proposed labels (hard stop).
#
# Run on Sherlock:  module load R/4.2.0; Rscript v2_G5_exvivo_B_subcluster.R
# Outputs (output/stepE/exvivo_B/, every file date-prefixed): UMAP + score + dot-plot PDFs, markers TSV,
#   cluster summary TSV with the proposed labels, per-cell table (cell, sample, condition, cluster, scores, proposed).
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

# Step 1. Files, folders and settings
F6_FILE        <- "${PROJECT_ROOT}/scRNA/2026-10-04_cHL_ex_vivo/output/2026-10-04_chl_exvivo_final_seurat_obj_f6.rds"
SIGNATURE_FILE <- "${PROJECT_ROOT}/cHL_v2/output/stepD/input/sigmtx_merged.txt"
OUTPUT_DIR     <- "${PROJECT_ROOT}/cHL_v2/output/stepE/exvivo_B"
DATE_PREFIX    <- "2026-10-25"
SEED           <- 42      # seed for PCA, Harmony, clustering and UMAP
N_HVG          <- 2000    # variable genes for the B-only sub-UMAP
N_PCS          <- 30      # principal components (author sub-UMAP recipe)
RESOLUTION     <- 0.4     # Louvain resolution (author sub-UMAP recipe)
N_SCORE_GENES  <- 30      # genes per preGC / postGC score
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_file <- function(name) file.path(OUTPUT_DIR, paste0(DATE_PREFIX, "_", name))

# Step 2. The B cells of the f6 object
exvivo <- readRDS(F6_FILE)
b_cells <- subset(exvivo, subset = cell_type == "B")
rm(exvivo)
invisible(gc())
DefaultAssay(b_cells) <- "RNA"
cat(sprintf("B cells: %d; per sample: %s\n", ncol(b_cells),
            paste(names(table(b_cells$sample)), table(b_cells$sample), sep = "=", collapse = " ")))

# Step 3. Sub-UMAP: LogNormalize -> HVG -> scale -> PCA -> Harmony on sample -> Louvain -> UMAP
set.seed(SEED)
b_cells <- NormalizeData(b_cells, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
b_cells <- FindVariableFeatures(b_cells, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
b_cells <- ScaleData(b_cells, verbose = FALSE)
b_cells <- RunPCA(b_cells, npcs = N_PCS, seed.use = SEED, verbose = FALSE)
b_cells <- RunHarmony(b_cells, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:N_PCS, verbose = FALSE)
b_cells <- FindNeighbors(b_cells, reduction = "harmony", dims = 1:N_PCS, verbose = FALSE)
b_cells <- FindClusters(b_cells, resolution = RESOLUTION, random.seed = SEED, verbose = FALSE)   # Louvain (default)
b_cells <- RunUMAP(b_cells, reduction = "harmony", dims = 1:N_PCS, seed.use = SEED, verbose = FALSE)

# Step 4. preGC / postGC genes from the round-14 signature (log2 ratio of signature values)
signature <- read.delim(SIGNATURE_FILE, row.names = 1, check.names = FALSE)
stopifnot(all(c("B_preGC", "B_GC", "B_Memory") %in% colnames(signature)))
postgc_mean <- (signature$B_GC + signature$B_Memory) / 2
log2_ratio  <- log2((signature$B_preGC + 1) / (postgc_mean + 1))
names(log2_ratio) <- rownames(signature)
log2_ratio  <- log2_ratio[names(log2_ratio) %in% rownames(b_cells)]          # genes present in the ex vivo panel
pregc_genes  <- names(sort(log2_ratio, decreasing = TRUE))[1:N_SCORE_GENES]
postgc_genes <- names(sort(log2_ratio, decreasing = FALSE))[1:N_SCORE_GENES]
cat("preGC genes :", paste(pregc_genes, collapse = " "), "\n")
cat("postGC genes:", paste(postgc_genes, collapse = " "), "\n")
b_cells <- AddModuleScore(b_cells, features = list(pregc_genes, postgc_genes), name = "score", seed = SEED)
b_cells$preGC_score  <- b_cells$score1
b_cells$postGC_score <- b_cells$score2

# Step 5. Markers by the author rule (rank by p_adj, then avg_log2FC; p_adj < 0.05; top 50 per cluster)
markers <- FindAllMarkers(b_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
markers <- as.data.table(markers)
setorder(markers, cluster, p_val_adj, -avg_log2FC)
top_markers <- markers[, head(.SD, 50), by = cluster]
top_markers[, passes := p_val_adj < 0.05]
fwrite(top_markers, out_file("exvivo_B_markers_top50.tsv"), sep = "\t")

# Step 6. Cluster summary and the PROPOSED label per cluster (higher mean score)
cluster_rows <- list()
for (this_cluster in levels(b_cells$seurat_clusters)) {      # not "cluster": that is a column of top_markers
  in_cluster <- b_cells$seurat_clusters == this_cluster
  mean_pregc  <- mean(b_cells$preGC_score[in_cluster])
  mean_postgc <- mean(b_cells$postGC_score[in_cluster])
  if (mean_pregc >= mean_postgc) {
    proposed <- "B_preGC"
  } else {
    proposed <- "B_postGC"
  }
  cluster_markers <- top_markers[as.character(top_markers$cluster) == this_cluster]
  cluster_rows[[length(cluster_rows) + 1]] <- data.table(
    cluster = this_cluster, cells = sum(in_cluster),
    per_sample = paste(names(table(b_cells$sample[in_cluster])), table(b_cells$sample[in_cluster]), sep = "=", collapse = " "),
    preGC_score = round(mean_pregc, 3), postGC_score = round(mean_postgc, 3), proposed = proposed,
    markers_passing_of_50 = sprintf("%d/%d", sum(cluster_markers$passes), nrow(cluster_markers)),
    top10_markers = paste(head(cluster_markers$gene, 10), collapse = ","))
}
cluster_summary <- rbindlist(cluster_rows)
fwrite(cluster_summary, out_file("exvivo_B_cluster_summary.tsv"), sep = "\t")
print(cluster_summary)

# Step 7. Per-cell table with the proposed label
b_cells$proposed_label <- cluster_summary$proposed[match(as.character(b_cells$seurat_clusters), cluster_summary$cluster)]
cell_table <- data.table(cell = colnames(b_cells), sample = b_cells$sample, condition = as.character(b_cells$condition),
                         cluster = as.character(b_cells$seurat_clusters), preGC_score = b_cells$preGC_score,
                         postGC_score = b_cells$postGC_score, proposed_label = b_cells$proposed_label)
fwrite(cell_table, out_file("exvivo_B_cells.tsv.gz"), sep = "\t")

# Step 8. Figures (PDF): UMAP by cluster / sample / proposed label; score feature plots; marker dot plot
umap_page <- DimPlot(b_cells, group.by = "seurat_clusters", label = TRUE) + ggtitle("Ex vivo B cells: clusters (res 0.4)") +
  DimPlot(b_cells, group.by = "sample") + ggtitle("Sample") +
  DimPlot(b_cells, group.by = "proposed_label") + ggtitle("Proposed label (higher score)")
ggsave(out_file("exvivo_B_umap.pdf"), umap_page, width = 16, height = 5)
score_page <- FeaturePlot(b_cells, features = c("preGC_score", "postGC_score"), ncol = 2) +
  VlnPlot(b_cells, features = c("preGC_score", "postGC_score"), group.by = "seurat_clusters", pt.size = 0, ncol = 2)
ggsave(out_file("exvivo_B_scores.pdf"), score_page, width = 12, height = 9)
known_markers <- c("IGHD", "TCL1A", "FCER2", "IL4R", "BACH2", "SELL",           # naive / pre-GC
                   "BCL6", "RGS13", "MEF2B", "AICDA", "LMO2", "S1PR2", "GCSAM",  # germinal centre
                   "CD27", "TNFRSF13B", "AIM2", "CD80", "SCIMP")                  # memory
known_markers <- known_markers[known_markers %in% rownames(b_cells)]
dot_page <- DotPlot(b_cells, features = known_markers, group.by = "seurat_clusters") + RotatedAxis() +
  ggtitle("Known naive/pre-GC, GC and memory markers by cluster")
ggsave(out_file("exvivo_B_marker_dotplot.pdf"), dot_page, width = 11, height = 4.5)
cat("outputs in", OUTPUT_DIR, "\n")
