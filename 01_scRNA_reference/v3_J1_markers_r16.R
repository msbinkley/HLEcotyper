# =============================================================================
# v3_J1_markers_r16.R — marker genes of the 16 round-16 cell types, to choose ONE representative marker panel for
# Fig S1c (author 2026-11-02: "why are there two fig S1C? you need to find the best representative markers and merge
# these into single figS1C").
#
# Author marker rule (CLAUDE.md): FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) on the
#   log-normalised RNA assay, no downsampling; rank by p_val_adj, then avg_log2FC; keep p_val_adj < 0.05; report the
#   top 50 per cell type and how many of the 50 pass.
# Dot-plot statistics for candidate genes, as v1 / v2_F1 (Seurat DotPlot): pct.exp = % of cells with expression > 0,
#   avg.exp = mean(expm1(log-normalised data)).
#
# Run on Sherlock (from cHL_v2/script/v3):  sbatch v3_J1_markers_r16.sh
# Inputs : OAK cHL_v2/output/cHL_merged_harmony_v2_r16.rds (= v3/data/input/scrna/cHL_merged_harmony_v2_r16.rds)
#          2026-10-26_round16_labels.tsv.gz (final_label per cell; = v3 p2("labels"))
# Outputs (OAK cHL_v2/output/r16/stepJ/):
#   2026-11-02_r16_FindAllMarkers_all.tsv        every marker FindAllMarkers returns
#   2026-11-02_r16_markers_top50.tsv             top 50 per type by the author rule, with the pass count
#   2026-11-02_r16_marker_candidates_dotstats.tsv  pct.exp / avg.exp per cell type for the candidate genes
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(future)
})

# Step 1. Files
V2_DIR <- "${PROJECT_ROOT}/cHL_v2"
OBJECT_FILE <- file.path(V2_DIR, "output/cHL_merged_harmony_v2_r16.rds")
LABEL_FILE <- file.path(V2_DIR, "output/r16/stepJ/2026-10-26_round16_labels.tsv.gz")   # copied from v3/data/input/scrna
OUT_DIR <- file.path(V2_DIR, "output/r16/stepJ")
SEED <- 20251102                                  # FindAllMarkers has no random step with presto; seeded anyway
TOP_N <- 50                                       # author rule: report the top 50 per cell type
CANDIDATES_PER_TYPE <- 30                         # genes per type carried into the dot-plot statistics

# Step 2. Object and final labels (one label per cell, 16 types)
object <- readRDS(OBJECT_FILE)
# fread() cannot open .gz on Sherlock (no R.utils; job 50845095 failed here), so read it with base R
labels <- read.delim(gzfile(LABEL_FILE), stringsAsFactors = FALSE)
labels <- as.data.table(labels)
print(head(labels, 3))
final_label <- setNames(labels$final_label, labels$cell)          # columns: cell, dataset, ..., final_label
stopifnot(setequal(names(final_label), colnames(object)))
object$final_label <- final_label[colnames(object)]
Idents(object) <- "final_label"
print(table(Idents(object)))

# Step 3. Log-normalised RNA data, one layer (the merged object may keep one layer per dataset)
DefaultAssay(object) <- "RNA"
object <- JoinLayers(object)
if (!"data" %in% Layers(object[["RNA"]])) {
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
}

# Step 4. FindAllMarkers with the author settings. Sherlock's R has no presto, so Seurat's own Wilcoxon test runs,
#         spread over 8 forked workers (same test and result, only faster)
plan("multicore", workers = 8)
options(future.globals.maxSize = 64 * 1024^3)    # 64 GB: the joined expression layer is passed to the workers
set.seed(SEED)
markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers <- as.data.table(markers)
fwrite(markers, file.path(OUT_DIR, "2026-11-02_r16_FindAllMarkers_all.tsv"), sep = "\t")

# Step 5. Top 50 per cell type: rank by p_val_adj, then by avg_log2FC (larger first); count how many pass p_adj < 0.05
top_rows <- list()
for (cell_type in sort(unique(as.character(markers$cluster)))) {
  one_type <- markers[cluster == cell_type]
  one_type <- one_type[order(p_val_adj, -avg_log2FC)]
  top <- head(one_type, TOP_N)
  top$rank <- seq_len(nrow(top))
  top$n_pass_in_top <- sum(top$p_val_adj < 0.05)
  top_rows[[cell_type]] <- top
  cat(sprintf("%-16s top %d: %d pass p_adj < 0.05\n", cell_type, nrow(top), sum(top$p_val_adj < 0.05)))
}
top_table <- rbindlist(top_rows)
fwrite(top_table, file.path(OUT_DIR, "2026-11-02_r16_markers_top50.tsv"), sep = "\t")

# Step 6. Candidate genes: per type, the passing markers with the largest specificity (pct.1 - pct.2), plus the
#         genes of the two existing S1c panels (v1 markers and the CP-A panel), so they can be compared
candidate_genes <- c()
for (cell_type in names(top_rows)) {
  passing <- markers[cluster == cell_type & p_val_adj < 0.05]
  passing <- passing[order(-(pct.1 - pct.2), -avg_log2FC)]
  candidate_genes <- c(candidate_genes, head(passing$gene, CANDIDATES_PER_TYPE))
}
existing_panels <- c("TNFRSF8", "CCL17", "CCL22", "IL13", "CD79A", "CD19", "MS4A1", "CR2", "FCER2", "AICDA", "GCSAM",
                     "ZBTB32", "TNFRSF13B", "XBP1", "SDC1", "CD3D", "CD4", "CCR7", "CD40LG", "FOXP3", "CTLA4",
                     "CD8A", "CXCR6", "FGFBP2", "GZMH", "CDCA5", "MKI67", "NKG7", "KLRF1", "XCL2", "IL7R", "RORC",
                     "FCGR1A", "MARCO", "CD68", "CD1E", "CLEC10A", "CLEC4C", "IL3RA", "CCN2", "HAS2", "COL1A1",
                     "COL3A1", "DCN", "CD274", "IL2RA", "TCL1A", "RGS13", "BCL6", "JCHAIN", "MZB1", "IRF4",
                     "PECAM1", "PTPRC", "CD3E", "CD27", "FCGR3A", "TOP2A", "GNLY", "KLRD1", "FOS", "JUN", "KIT",
                     "CD163", "CD14", "C1QA", "LYZ", "CD1C", "LILRA4", "COL1A2", "EGFL7", "VCAN", "HSPA1A",
                     "PDGFRA", "VWF", "TRDC", "TRGC1", "LAMP3", "CD8B", "GZMK", "FCRL4", "PRDM1")
candidate_genes <- unique(c(candidate_genes, existing_panels))
candidate_genes <- candidate_genes[candidate_genes %in% rownames(object)]
cat("candidate genes:", length(candidate_genes), "\n")

# Step 7. Dot-plot statistics per cell type for the candidates (as Seurat DotPlot / v2_F1)
expression <- LayerData(object, assay = "RNA", layer = "data")[candidate_genes, , drop = FALSE]
cell_types <- sort(unique(as.character(Idents(object))))
dot_rows <- list()
for (cell_type in cell_types) {
  cells <- which(as.character(Idents(object)) == cell_type)
  one_type <- expression[, cells, drop = FALSE]
  pct_exp <- Matrix::rowMeans(one_type > 0) * 100
  avg_exp <- Matrix::rowMeans(expm1(one_type))
  dot_rows[[cell_type]] <- data.table(group = cell_type, feature = candidate_genes, pct.exp = pct_exp,
                                      avg.exp = avg_exp, n_cells = length(cells))
}
fwrite(rbindlist(dot_rows), file.path(OUT_DIR, "2026-11-02_r16_marker_candidates_dotstats.tsv"), sep = "\t")
cat("written to", OUT_DIR, "\n")
