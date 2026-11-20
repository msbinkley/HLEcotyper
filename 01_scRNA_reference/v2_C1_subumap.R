# =============================================================================
# v2_C1_subumap.R — one lineage sub-UMAP of round r2 (Sherlock); the major type is given by env var TYPE.
#
# Recipe = author's ex vivo recluster():
#   SCTransform v2 (glmGamPoi) refit on the type's cells -> PCA 30 -> Harmony (dataset + chemistry + sample, seed 42)
#   -> Louvain res 0.4 (Seurat default, as in the ex vivo recluster()) -> UMAP.
#   Markers on log-normalised RNA (Wilcoxon).
# Input : output/cHL_merged_v2_r2_qc.rds (from v2_C0_qc_scores.R)
# Output (output/subumap/r2_sub_<TYPE>_*):
#   obj.rds                clustered object (checkpoint: reused on a rerun)
#   cell_ids.rds           cell id + subcluster (read by v2_C2_recluster_r3.R)
#   cells.tsv.gz           cell ids + UMAP coordinates + QC flags
#   top50_markers.tsv      top-50 markers per subcluster
#   qc_per_subcluster.tsv  per-subcluster QC table
#   umaps.pdf              UMAP PDFs
# Run   : TYPE=<major type> Rscript v2_C1_subumap.R      (e.g. TYPE=Myeloid; one job per major type)
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
  library(ggplot2)
  library(patchwork)
})
TYPE <- Sys.getenv("TYPE")
stopifnot(nzchar(TYPE))                                             # TYPE must be set
options(future.globals.maxSize = 16 * 1024^3)   # SCTransform exports >1 GB of globals for the large types
W <- "${PROJECT_ROOT}/cHL_v2"             # project folder on Sherlock
OUT <- file.path(W, "output/subumap")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))
tag <- paste0("r2_sub_", TYPE)

# out_file(): output path for this type, e.g. out_file("obj.rds") -> output/subumap/r2_sub_<TYPE>_obj.rds
out_file <- function(name) file.path(OUT, paste0(tag, "_", name))

# log_(): print a time-stamped message
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")


# ---- Step 1. Sub-UMAP of this major type (or reload it if it was already made) -----------------------------
if (file.exists(out_file("obj.rds"))) {                              # checkpoint: reuse the clustered object on a rerun
  sub <- readRDS(out_file("obj.rds"))
  log_("%s: reloaded clustered object (%d cells)", TYPE, ncol(sub))
} else {
  # 1a. the cells of this major type, without the whole-object reductions / graphs
  obj <- readRDS(file.path(W, "output/cHL_merged_v2_r2_qc.rds"))
  sub <- subset(obj, cells = colnames(obj)[obj$major_type == TYPE])
  rm(obj)
  invisible(gc())
  sub@reductions <- list()
  sub@graphs <- list()
  sub@neighbors <- list()
  log_("%s: %d cells", TYPE, ncol(sub))

  # 1b. SCTransform v2 -> PCA 30 -> unintegrated UMAP (seed 42) -> Harmony (seed 42) -> clusters res 0.4 (seed 42) -> UMAP (seed 42)
  sub <- SCTransform(sub, vst.flavor = "v2", method = "glmGamPoi", verbose = FALSE)
  sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
  sub <- RunUMAP(sub, reduction = "pca", dims = 1:30, reduction.name = "umap_pca", reduction.key = "umappca_",
                 seed.use = 42, verbose = FALSE)
  set.seed(42)                                                      # seed for Harmony
  sub <- RunHarmony(sub, group.by.vars = c("dataset", "chemistry", "sample"), reduction.use = "pca",
                    reduction.save = "harmony", verbose = FALSE)
  sub <- FindNeighbors(sub, reduction = "harmony", dims = 1:30, verbose = FALSE)
  sub <- FindClusters(sub, resolution = 0.4, random.seed = 42, verbose = FALSE)
  sub <- RunUMAP(sub, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)
  sub$subcluster <- as.character(sub$seurat_clusters)
  saveRDS(sub, out_file("obj.rds"))
}
log_("%d subclusters", length(unique(sub$subcluster)))


# ---- Step 2. Markers (log-normalised RNA, Wilcoxon) -----------------------------------------------------------
# positive markers only; log2FC >= 0.25; expressed in >= 10% of cells; at most 500 cells per subcluster (seed 42)
DefaultAssay(sub) <- "RNA"
Idents(sub) <- "subcluster"
markers <- FindAllMarkers(sub, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1, max.cells.per.ident = 500,
                          random.seed = 42, verbose = FALSE)
markers <- data.table(markers)
markers[, symbol := v47$symbol[match(gene, v47$ensg)]]

# number of significant markers (adjusted p < 0.05) per subcluster
n_significant <- markers[p_val_adj < 0.05, .(n_sig_markers = .N), by = .(subcluster = as.character(cluster))]

# top 50 significant markers per subcluster, ranked by log2 fold change
significant <- markers[p_val_adj < 0.05]
significant <- significant[order(cluster, -avg_log2FC)]
top <- significant[, head(.SD, 50), by = cluster]
fwrite(top[, .(subcluster = cluster, symbol, avg_log2FC = round(avg_log2FC, 3), pct.1, pct.2, p_val_adj)],
       out_file("top50_markers.tsv"), sep = "\t")


# ---- Step 3. QC per subcluster -----------------------------------------------------------------------------------
md <- data.table(sub@meta.data, keep.rownames = "cell")             # per-cell table

# 3a. the 3 most common published labels per subcluster with their %, e.g. "T_CD4 80; T_Treg 15; S1prov 5"
#     ("S1prov" = S1 cells with only a provisional label)
label_counts <- md[, .N, by = .(subcluster, label)]
label_counts <- label_counts[order(-N)]
label_counts[, label_name := fifelse(is.na(label), "S1prov", label)]
label_counts[, label_pct := round(100 * N / sum(N)), by = subcluster]
lab_comp <- label_counts[, .(top_labels = paste(head(paste0(label_name, " ", label_pct), 3), collapse = "; ")), by = subcluster]

# 3b. cells per dataset per subcluster (one column per dataset)
dataset_counts <- dcast(md[, .N, by = .(subcluster, dataset)], subcluster ~ dataset, value.var = "N", fill = 0)

# 3c. top_contam(): the 1-2 most frequently flagged foreign lineages among a subcluster's cells ("" = none)
top_contam <- function(contam_lineages) {
  flagged <- contam_lineages[contam_lineages != ""]
  lineage_counts <- table(unlist(strsplit(flagged, ",")))
  if (length(lineage_counts)) {
    sorted_names <- names(sort(lineage_counts, decreasing = TRUE))
    paste(sorted_names[1:min(2, length(lineage_counts))], collapse = ",")
  } else {
    ""
  }
}

# 3d. per subcluster: cells, samples, median UMI / genes / % MT / % ribo (original counts), % flagged as contaminated,
#     protected T_gd / HRS cells, and the most common flagged lineages
qc_table <- md[, .(cells = .N,
                   n_samples = uniqueN(sample),
                   med_nCount_orig = as.numeric(median(nCount_orig)),
                   med_nFeature_orig = as.numeric(median(nFeature_orig)),
                   med_pct_mt = round(median(percent.mt), 2),
                   med_pct_ribo = round(median(percent.ribo), 2),
                   pct_contaminated = round(100 * mean(contaminated), 1),
                   n_protected_Tgd = sum(label %in% "T_gd"),             # rare: caution before removal
                   n_protected_HRS = sum(label %in% "HRS"),              # rare: caution before removal
                   top_contam = top_contam(contam_lineages)),
               by = subcluster]

# 3e. median of every lineage / stress score per subcluster (3 decimals)
score_cols <- grep("_score$", names(md), value = TRUE)
score_medians <- md[, lapply(.SD, function(v) round(median(v), 3)), by = subcluster, .SDcols = score_cols]

# 3f. join everything by subcluster, sort by subcluster number, save
qc_table <- merge(qc_table, score_medians, by = "subcluster")
qc_table <- merge(qc_table, dataset_counts, by = "subcluster", all.x = TRUE)
qc_table <- merge(qc_table, lab_comp, by = "subcluster", all.x = TRUE)
qc_table <- merge(qc_table, n_significant, by = "subcluster", all.x = TRUE)
qc_table <- qc_table[order(as.integer(subcluster))]
fwrite(qc_table, out_file("qc_per_subcluster.tsv"), sep = "\t")

# 3g. per-cell table with the Harmony UMAP coordinates
umap_harmony <- Embeddings(sub, "umap")
umap_pca <- Embeddings(sub, "umap_pca")                             # (not used below; kept from the original)
cell_export <- md[, .(cell, dataset, sample, label, r2_cluster, subcluster, contaminated, contam_lineages, protected)]
cell_export[, UMAP_1 := umap_harmony[, 1]]
cell_export[, UMAP_2 := umap_harmony[, 2]]
fwrite(cell_export, out_file("cells.tsv.gz"), sep = "\t")


# ---- Step 4. UMAP PDFs ----------------------------------------------------------------------------------------------
# label shown on the plots: published label, or "S1prov:<type>" for S1 cells with only a provisional label
lab_show <- ifelse(is.na(sub$label), paste0("S1prov:", sub$s1_cell_type), sub$label)
sub$label_show <- lab_show
## reduction pinned: Seurat's default lookup matches "umap_pca" (computed first) before the Harmony "umap"
pdf(out_file("umaps.pdf"), width = 11, height = 9)

# page 1: subclusters
title_text <- sprintf("%s: %d cells, %d subclusters (res 0.4)", TYPE, ncol(sub), length(unique(sub$subcluster)))
print(DimPlot(sub, reduction = "umap", group.by = "subcluster", label = TRUE, repel = TRUE, pt.size = 0.2) +
        ggtitle(title_text))
# page 2: published labels
print(DimPlot(sub, reduction = "umap", group.by = "label_show", label = TRUE, repel = TRUE, pt.size = 0.2) +
        ggtitle(paste(TYPE, ": published labels")))
# page 3: dataset
print(DimPlot(sub, reduction = "umap", group.by = "dataset", shuffle = TRUE, pt.size = 0.2) +
        ggtitle(paste(TYPE, ": dataset")))
# page 4: subclusters split by dataset
print(DimPlot(sub, reduction = "umap", group.by = "subcluster", split.by = "dataset", ncol = 3, pt.size = 0.2) +
        NoLegend() +
        ggtitle(paste(TYPE, ": subclusters split by dataset")))
# page 5: cells flagged as contaminated (red)
print(DimPlot(sub, reduction = "umap", group.by = "contaminated", shuffle = TRUE, pt.size = 0.2,
              cols = c(`FALSE` = "grey85", `TRUE` = "red")) +
        ggtitle(paste(TYPE, ": cells with a high foreign-lineage score")))
# page 6: lineage / stress scores
print(FeaturePlot(sub, reduction = "umap", features = score_cols, ncol = 3, pt.size = 0.1, order = TRUE) &
        theme(plot.title = element_text(size = 9)))
# page 7: QC violins per subcluster
print(VlnPlot(sub, features = c("nCount_orig", "nFeature_orig", "percent.mt", "Stress_score"), group.by = "subcluster",
              pt.size = 0, ncol = 2))
# page 8: canonical markers per subcluster (genes shown by symbol)
dot_symbols <- c("CD3E", "CD4", "CD8A", "FOXP3", "MKI67", "NKG7", "GNLY", "TRDC", "KIT", "MS4A1", "CD79A", "BCL6", "CD27", "MZB1",
                 "JCHAIN", "TNFRSF8", "CCL17", "LYZ", "CD68", "CD1C", "CLEC9A", "LILRA4", "COL1A2", "HSPA1A", "FOS", "MT1X", "MT-ND2")
dot_ensg <- v47$ensg[match(dot_symbols, v47$symbol)]
keep <- !is.na(dot_ensg) & dot_ensg %in% rownames(sub)              # only genes that are in the object
print(DotPlot(sub, assay = "RNA", features = dot_ensg[keep], group.by = "subcluster") +
        scale_x_discrete(labels = setNames(dot_symbols[keep], dot_ensg[keep])) +
        RotatedAxis() +
        ggtitle(paste(TYPE, ": canonical markers per subcluster")))
# page 9: unintegrated (PCA) UMAP by dataset
print(DimPlot(sub, reduction = "umap_pca", group.by = "dataset", shuffle = TRUE, pt.size = 0.2) +
        ggtitle(paste(TYPE, ": unintegrated (PCA) UMAP")))
dev.off()


# ---- Step 5. Cell ids + subcluster (input of v2_C2_recluster_r3.R) ----------------------------------------------------
saveRDS(data.table(cell = colnames(sub), subcluster = sub$subcluster), out_file("cell_ids.rds"))
log_("done %s", TYPE)
