# =============================================================================
# v2_C7_subumap_r12_cluster.R — sub-UMAP of one round-12 cluster (Sherlock, or locally with the env overrides below).
#
# PURPOSE
#   env CLUSTER (default 8 = T_MKI67: cycling T + cycling B_GC; author 2026-10-20 "perform sub-UMAP analysis on 8 T_MKI67
#   and annotate its subtypes").
#   Recipe = cluster 4/8/16 sub-UMAP (v2_C5f, author-approved):
#     LogNormalize -> 3,000 variable genes per dataset (ranked by # datasets; datasets with < 50 cells skipped for HVG)
#     -> scale -> PCA 30 -> Harmony (dataset + chemistry + sample, theta 4/4/2, seed 42) -> Louvain res 0.4 / 0.8 / 1.2
#     -> UMAP (seed 42).
#   Markers (author rule): FindAllMarkers on log-normalised RNA, only.pos, min.pct 0.25, logfc 0.25, no downsampling;
#     top 50 per subcluster ranked by p_val_adj then avg_log2FC; how many of the 50 have p_val_adj < 0.05.
#     Resolutions 0.8 and 1.2.
#   QC per subcluster: cells, samples, datasets, r12 / published labels, median UMI vs dataset median, % mt, contamination,
#     stress and lineage scores (cpB), % cells expressing key genes.
#
# ENV VARIABLES (paths can be overridden so the same script runs locally, author 2026-10-20)
#   W       = cHL_v2 root on OAK                    OUT    = output folder (default W/output/subumap)
#   OBJ     = round-12 object                       LABELS = round-12 labels file
#   QC      = per-cell QC table (cpB)               GENES  = v47 gene table
#   CLUSTER = round-12 res-0.6 cluster to take      ROUND  = "r2", "r3", ... for cleanup rounds (added to file names)
#   Cleanup round (author 2026-10-20): drop subclusters DROP (comma list) of column DROP_COL in the previous round's
#   cells file DROP_FROM; N_EXPECT = expected number of dropped cells.
#     cl8 round 2: DROP_FROM=<r12_cl8_cells.tsv.gz> DROP_COL=sub_res0.8 DROP=10,12 N_EXPECT=275
#                  (low UMI: 0.36 / 0.53 of dataset median)
#
# OUTPUT : OUT/r12_cl<CLUSTER>[_<ROUND>]_* (object, cells, markers, QC tables, PDFs; removed_cells in cleanup rounds)
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(data.table)
  library(harmony)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
options(future.globals.maxSize = 16 * 1024^3)

W      <- Sys.getenv("W", "${PROJECT_ROOT}/cHL_v2")
OUT    <- Sys.getenv("OUT", file.path(W, "output/subumap"))
OBJ    <- Sys.getenv("OBJ", file.path(W, "output/cHL_merged_harmony_v2_r12.rds"))
LABELS <- Sys.getenv("LABELS", file.path(W, "output/stepC/inputs/2026-10-20_round12_labels.tsv.gz"))
QC     <- Sys.getenv("QC", file.path(W, "output/cpB_cell_qc.tsv.gz"))
GENES  <- Sys.getenv("GENES", file.path(W, "ref/gene_table_v47.tsv"))
CL     <- Sys.getenv("CLUSTER", "8")
ROUND  <- Sys.getenv("ROUND", "")       # ROUND = "r2", "r3", ... for cleanup rounds

# Output file names start with r12_cl<CLUSTER>_ (plus <ROUND>_ in cleanup rounds)
if (nzchar(ROUND)) {
  tag <- paste0("r12_cl", CL, "_", ROUND)
} else {
  tag <- paste0("r12_cl", CL)
}
out_path <- function(x) file.path(OUT, paste0(tag, "_", x))

DROP_FROM <- Sys.getenv("DROP_FROM", "")
DROP_COL  <- Sys.getenv("DROP_COL", "sub_res0.8")
DROP      <- strsplit(Sys.getenv("DROP", ""), ",")[[1]]
N_EXPECT  <- as.integer(Sys.getenv("N_EXPECT", "0"))
v47 <- fread(GENES)

# Print a time-stamped log line
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")


## Step 1. Cells of the chosen cluster (minus dropped subclusters in a cleanup round) -------------------------------------------
full_obj   <- readRDS(OBJ)
r12_labels <- fread(cmd = sprintf("gzip -dc '%s'", LABELS))
keep <- colnames(full_obj)[as.character(full_obj$RNA_snn_res.0.6) == CL]
if (nzchar(DROP_FROM)) {
  previous_round <- fread(cmd = sprintf("gzip -dc '%s'", DROP_FROM))
  drop <- previous_round[as.character(get(DROP_COL)) %in% DROP, cell]      # cells in the subclusters to drop
  stopifnot(length(drop) == N_EXPECT, all(drop %in% keep))
  keep <- setdiff(keep, drop)
  removed_subcluster <- previous_round[match(drop, cell), get(DROP_COL)]
  fwrite(data.table(cell = drop, removed_subcluster = removed_subcluster), out_path("removed_cells.tsv.gz"), sep = "\t")
  log_("dropped %d cells (%s %s)", length(drop), DROP_COL, paste(DROP, collapse = ","))
}

# Step 2. New Seurat object: raw counts (rows renamed to unique v47 symbols) + metadata + round-12 label
cnt <- LayerData(full_obj, assay = "RNA", layer = "counts")[, keep]
md <- full_obj@meta.data[keep, c("dataset", "sample", "chemistry", "label")]
rm(full_obj)
invisible(gc())
md$r12_label <- r12_labels$r12_label[match(keep, r12_labels$cell)]
rownames(cnt) <- v47$symbol[match(rownames(cnt), v47$ensg)]
stopifnot(!anyNA(rownames(cnt)), !anyDuplicated(rownames(cnt)))
sub <- CreateSeuratObject(cnt, meta.data = md)
rm(cnt)
invisible(gc())
log_("cluster %s: %d cells; %s", CL, ncol(sub), paste(names(table(sub$dataset)), table(sub$dataset), sep = "=", collapse = " "))


## Step 3. Normalise, variable genes, PCA, Harmony, clustering, UMAP ----------------------------------------------------------
sub <- NormalizeData(sub, verbose = FALSE)

# 3a. 3,000 variable genes per dataset, only for datasets with >= 50 cells in this cluster
cells_per_dataset <- table(sub$dataset)
use <- names(cells_per_dataset)[cells_per_dataset >= 50]
cells_by_dataset <- split(colnames(sub), sub$dataset)[use]
hv <- list()                                                     # one HVG vector per dataset
for (ds in names(cells_by_dataset)) {
  ds_obj <- CreateSeuratObject(LayerData(sub, assay = "RNA", layer = "counts")[, cells_by_dataset[[ds]]])
  ds_obj <- FindVariableFeatures(ds_obj, nfeatures = 3000, verbose = FALSE)
  hv[[ds]] <- VariableFeatures(ds_obj)
}
rm(ds_obj)

# 3b. Rank genes by how many datasets pick them, then by their median rank within those datasets; keep the top 3,000
rk <- table(unlist(hv))                                          # number of datasets picking each gene
rank_matrix <- matrix(NA_integer_, nrow = length(rk), ncol = length(hv))   # gene x dataset: rank in that dataset's HVG list
for (d in seq_along(hv)) {
  rank_matrix[, d] <- match(names(rk), hv[[d]])
}
med <- apply(rank_matrix, 1, median, na.rm = TRUE)               # median rank over the datasets that pick the gene
hvg <- names(rk)[order(-rk, med)][1:3000]
VariableFeatures(sub) <- hvg

# 3c. Scale, PCA 30 (seed 42), Harmony theta 4/4/2 (seed 42), Louvain res 0.4 / 0.8 / 1.2 (seed 42), UMAP (seed 42)
sub <- ScaleData(sub, features = hvg, verbose = FALSE)
sub <- RunPCA(sub, features = hvg, npcs = 30, seed.use = 42, verbose = FALSE)
set.seed(42)                                                     # Harmony uses R's random numbers
sub <- RunHarmony(sub, group.by.vars = c("dataset", "chemistry", "sample"), theta = c(4, 4, 2),
                  reduction.use = "pca", reduction.save = "harmony", verbose = FALSE)
sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (r in c(0.4, 0.8, 1.2)) {
  sub <- FindClusters(sub, resolution = r, random.seed = 42, verbose = FALSE)
}
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)
saveRDS(sub, out_path("obj.rds"))
fwrite(data.table(cell = colnames(sub),
                  sub_res0.4 = sub$RNA_snn_res.0.4,
                  sub_res0.8 = sub$RNA_snn_res.0.8,
                  sub_res1.2 = sub$RNA_snn_res.1.2,
                  Embeddings(sub, "umap")),
       out_path("cells.tsv.gz"), sep = "\t")


## Step 4. QC + markers per resolution ---------------------------------------------------------------------------------------
qc <- fread(cmd = sprintf("gzip -dc '%s'", QC), na.strings = c("", "NA"))
qc <- qc[cell %in% colnames(sub)]
dataset_median_umi <- qc[, .(dsUMI = as.numeric(median(nCount_orig))), by = dataset]     # median UMI of each dataset (these cells)

# key genes for the "% cells expressing" columns (those present in the object)
G <- c("MKI67", "TOP2A", "CD3E", "CD4", "CD8A", "FOXP3", "IL2RA", "CXCL13", "PDCD1", "NKG7", "GZMB", "MS4A1", "CD19", "BCL6",
       "AICDA", "RGS13", "MZB1", "JCHAIN", "TNFRSF8", "IRF4")
G <- G[G %in% rownames(sub)]
expr <- as.matrix(LayerData(sub, assay = "RNA", layer = "counts")[G, ] > 0)          # gene x cell: expressed (count > 0)

# Helper: the two most common values and their % (e.g. "T_CD4 60; B_GC 30")
top2 <- function(x) {
  counts <- sort(table(x), decreasing = TRUE)
  n_show <- min(2, length(counts))
  paste(names(counts)[1:n_show], round(100 * counts[1:n_show] / length(x)), collapse = "; ")
}
# Helper: % of TRUE values, rounded
pct_true <- function(v) round(100 * mean(v))

for (r in c("0.8", "1.2")) {
  Idents(sub) <- paste0("RNA_snn_res.", r)

  # 4a. Markers (author rule): top 50 per subcluster by p_val_adj, then avg_log2FC
  markers_all <- FindAllMarkers(sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  top50 <- markers_all %>%
    group_by(cluster) %>%
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 50)
  fwrite(top50, out_path(sprintf("top50_markers_res%s.tsv", r)), sep = "\t")
  # per subcluster: how many of the 50 have p_adj < 0.05, the largest p_adj, and the first 15 genes
  marker_summary <- data.table(top50)[, .(n_sig    = sum(p_val_adj < 0.05),
                                          max_padj = signif(max(p_val_adj), 2),
                                          top15    = paste(head(gene, 15), collapse = " ")),
                                      by = .(subcluster = as.character(cluster))]

  # 4b. Per-cell table: subcluster, labels, QC values, UMI relative to the dataset median
  cell_qc <- data.table(cell = colnames(sub), subcluster = as.character(Idents(sub)), dataset = sub$dataset,
                        sample = sub$sample, label = sub$label, r12_label = sub$r12_label)
  cell_qc <- merge(cell_qc, qc[, .(cell, nCount_orig, percent.mt, contaminated, Stress_score, T_score, B_score, Plasma_score, HRS_score)],
                   by = "cell")
  cell_qc <- merge(cell_qc, dataset_median_umi, by = "dataset")
  cell_qc[, relUMI := nCount_orig / dsUMI]

  # 4c. % of each subcluster's cells expressing each key gene (columns pct_<gene>)
  expr_table <- data.table(cell = colnames(sub), t(expr))
  expr_table <- expr_table[cell_qc[, .(cell, subcluster)], on = "cell"]          # rows in cell_qc order, + subcluster
  pct_expressing <- expr_table[, lapply(.SD, pct_true), by = subcluster, .SDcols = G]
  setnames(pct_expressing, G, paste0("pct_", G))

  # 4d. One QC row per subcluster (medians; published label: "S1" for S1 cells)
  qc_table <- cell_qc[, .(cells      = .N,
                          samples    = uniqueN(sample),
                          relUMI     = round(median(relUMI), 2),
                          mt         = round(median(percent.mt), 1),
                          contam     = round(100 * mean(contaminated), 1),
                          stress     = round(median(Stress_score), 2),
                          T          = round(median(T_score), 2),
                          B          = round(median(B_score), 2),
                          Plasma     = round(median(Plasma_score), 2),
                          HRS        = round(median(HRS_score), 2),
                          r12_labels = top2(r12_label),
                          published  = top2(ifelse(dataset == "s1", "S1", label)),
                          datasets   = paste(names(table(dataset)), table(dataset), sep = "=", collapse = " ")),
                      by = subcluster]
  qc_table <- merge(qc_table, marker_summary, by = "subcluster")
  qc_table <- merge(qc_table, pct_expressing, by = "subcluster")
  qc_table <- qc_table[order(as.integer(subcluster))]
  fwrite(qc_table, out_path(sprintf("qc_res%s.tsv", r)), sep = "\t")
  cat("\n== res", r, "\n")
  print(qc_table[, .(subcluster, cells, samples, relUMI, n_sig, r12_labels, pct_MKI67, pct_CD3E, pct_MS4A1, top15)])
}


## Step 5. Figures ---------------------------------------------------------------------------------------------------------------
cols <- c(B_GC = "#E6AB02", T_MKI67 = "#FADD4B", Plasma = "#FB9A99", HRS = "#1B9E77", T_CD4 = "#B2DF8A", T_Treg = "#A6CEE3", T_CD8 = "#999999",
          B_preGC = "#1F78B4", B_Memory = "#FF7F00", S1 = "#5AAE61", NK = "#FDBF6F", T_gd = "#A65628", ILC = "#00C19A")
sub$published <- ifelse(sub$dataset == "s1", "S1", sub$label)          # S1 cells have no published label
pdf(out_path("umaps.pdf"), width = 20, height = 8)
# page 1: published label | subclusters res 0.8 | subclusters res 1.2
print(DimPlot(sub, group.by = "published", shuffle = TRUE, pt.size = 0.8, cols = cols, label = TRUE, repel = TRUE) +
        ggtitle(sprintf("r12 cluster %s sub-UMAP (%d cells): published label", CL, ncol(sub))) |
      DimPlot(sub, group.by = "RNA_snn_res.0.8", pt.size = 0.8, label = TRUE) + ggtitle("subclusters res 0.8") |
      DimPlot(sub, group.by = "RNA_snn_res.1.2", pt.size = 0.8, label = TRUE) + ggtitle("subclusters res 1.2"))
# page 2: dataset | res 0.4
print(DimPlot(sub, group.by = "dataset", shuffle = TRUE, pt.size = 0.8) + ggtitle("dataset") |
      DimPlot(sub, group.by = "RNA_snn_res.0.4", pt.size = 0.8, label = TRUE) + ggtitle("res 0.4"))
# page 3: key genes on the UMAP
print(FeaturePlot(sub, features = G, ncol = 5, pt.size = 0.4, order = TRUE) &
        theme(plot.title = element_text(size = 9), axis.text = element_blank()))
# page 4: dot plots of the key genes at res 0.8 and 1.2
print(DotPlot(sub, features = G, group.by = "RNA_snn_res.0.8") + RotatedAxis() + ggtitle("res 0.8") |
      DotPlot(sub, features = G, group.by = "RNA_snn_res.1.2") + RotatedAxis() + ggtitle("res 1.2"))
dev.off()
log_("done")
