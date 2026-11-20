# =============================================================================
# v2_B_merge_harmony.R — step B of REANALYSIS_PLAN_v2: merge the five scRNA datasets and integrate (round 1).
#
# Runs on Sherlock (R/4.2.0, Seurat 5.5 + harmony 2.0.5).
#
# What it does:
#   - merges the five scRNA datasets on the step-A gene universe (raw counts, v47 ENSG rows),
#   - keeps the published labels,
#   - integrates with Harmony,
#   - exports slim QC tables for the author's UMAP review (checkpoint CP-A).
#
#   dataset  object                               label column (figure)   notes
#   lab      data/01_input/cHL_scRNA_seurat.rds   cell_annotation (Fig 2A) sample: neutral codes; chemistry 3prime / Flex_v1.0
#   aoki     2025-08-17_Aoki_scRNA_filtered.rds   final_2025_08_12 (Fig S2B) sample: StudyID (public)
#   stewart  external_Stewart_filtered.rds        final_2025_08_12 (Fig S2C) raw counts rebuilt (B0)
#   tenx     external_10X_filtered.rds            final_cluster2 (Fig S2A)   B_cell -> B_preGC, flagged ref_exclude
#   s1       ex vivo f6 object, Sample1           cell_type mapped; Monocyte / Innate_T / Endothelial / IEG_stress kept
#                                                 with label NA + provisional label (resolved at CP-A)
# Nothing is excluded from the object. ref_exclude only marks cells that step C will leave out of the reference.
#
# Recipe: LogNormalize -> 3,000 HVGs chosen per dataset and ranked by how many datasets pick them -> ScaleData ->
#         PCA 50 -> Harmony (dataset + chemistry) on 50 PCs -> neighbours on 30 Harmony dims ->
#         Louvain clusters at res 0.3 / 0.6 / 1.0 -> UMAP. Seeds: 42 (PCA, UMAP, Harmony), 0 (clustering).
#
# Inputs : <W>/input/*.rds, <W>/input/lab_cell_sample_map.tsv, the ex vivo f6 object (S1_OBJ),
#          <W>/ref/gene_universe.tsv, <W>/ref/gene_table_v47.tsv, <W>/ref/feature_map_<dataset>.tsv (step A)
# Outputs (in <W>/output):
#   cHL_merged_harmony_v2.rds                 merged + integrated object
#   hvg_3000.tsv                              the 3,000 HVGs
#   cpA_cells.tsv.gz                          per-cell table (metadata, clusters, UMAP coordinates, mixing scores)
#   cpA_dotplot_by_cluster_res0.6.tsv         marker mean expression / % expressing per res-0.6 cluster
#   cpA_dotplot_by_label.tsv                  the same per published label
#   cpA_sessionInfo.txt
# Run: Rscript v2_B_merge_harmony.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(harmony)
})
set.seed(42)                                                          # global seed

# ---- Settings ------------------------------------------------------------------------------------------
W   <- "${PROJECT_ROOT}/cHL_v2"           # project folder on Sherlock
REF <- file.path(W, "ref")
OUT <- file.path(W, "output")
dir.create(OUT, showWarnings = FALSE)
S1_OBJ <- "${PROJECT_ROOT}/scRNA/2026-10-04_cHL_ex_vivo/output/2026-10-04_chl_exvivo_final_seurat_obj_f6.rds"
universe <- fread(file.path(REF, "gene_universe.tsv"))$ensg
v47 <- fread(file.path(REF, "gene_table_v47.tsv"))

# log_(): print a time-stamped message
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

# to_universe(): features x cells counts -> universe genes x cells counts.
#   Each feature is mapped to its v47 gene with the step-A feature map; features that map to the same
#   universe gene are summed; features outside the universe are dropped.
to_universe <- function(counts, dataset) {
  feature_map <- fread(file.path(REF, sprintf("feature_map_%s.tsv", dataset)))
  feature_map <- feature_map[match(rownames(counts), feature_map$feature)]       # one map row per count row
  in_universe <- !is.na(feature_map$v47) & feature_map$v47 %in% universe
  # features x universe indicator matrix: 1 where a feature maps to a universe gene
  mapping <- sparseMatrix(i = which(in_universe), j = match(feature_map$v47[in_universe], universe), x = 1,
                          dims = c(nrow(counts), length(universe)))
  universe_counts <- crossprod(mapping, counts)                                   # = t(mapping) %*% counts
  rownames(universe_counts) <- universe
  log_("%-8s %6d cells: %5d features -> %5d universe genes; %.1f%% of UMIs kept", dataset, ncol(counts),
       nrow(counts), sum(rowSums(mapping) > 0 & TRUE), 100 * sum(universe_counts) / sum(counts))
  as(universe_counts, "CsparseMatrix")
}

parts <- list()          # counts per dataset (universe genes x cells)
meta  <- list()          # per-cell metadata per dataset


# ---- Step 1. lab ------------------------------------------------------------------------------------------
lab_obj <- readRDS(file.path(W, "input/lab_scRNA_seurat.rds"))
counts <- LayerData(lab_obj, assay = "RNA", layer = "counts")
cell_meta <- lab_obj@meta.data
rm(lab_obj)
invisible(gc())

# sample of each lab cell: from the barcode prefix (e.g. "scRNA_1"), via lab_cell_sample_map.tsv
sample_map <- fread(file.path(W, "input/lab_cell_sample_map.tsv"))               # its ID column uses "." for "-"
cell_ids <- gsub("-", ".", colnames(counts))
barcode_key <- sub("^(.+_)[ACGTN]+[.-]", "\\1", cell_ids)                         # e.g. "scRNA_1"
key_to_sample <- unique(data.table(bc_key = sub("^(.+_)[ACGTN]+[.-]", "\\1", sample_map$ID), sample = sample_map$Sample))
stopifnot(!anyDuplicated(key_to_sample$bc_key))                                  # barcode suffix -> sample is 1:1
lab_sample <- key_to_sample$sample[match(barcode_key, key_to_sample$bc_key)]
stopifnot(!anyNA(lab_sample))

parts$lab <- to_universe(counts, "lab")
meta$lab <- data.table(cell = paste0("lab_", colnames(counts)), dataset = "lab", sample = lab_sample,
                       chemistry = ifelse(cell_meta$orig.ident == "snRNA", "Flex_v1.0", "3prime"),
                       label = sub("^T_MKI67\\+$", "T_MKI67", cell_meta$cell_annotation),     # unify "T_MKI67+"
                       label_source = "published (Fig 2A)")


# ---- Step 2. Aoki -------------------------------------------------------------------------------------------
aoki_obj <- readRDS(file.path(W, "input/aoki_filtered.rds"))
counts <- LayerData(aoki_obj, assay = "originalexp", layer = "counts")
cell_meta <- aoki_obj@meta.data
rm(aoki_obj)
invisible(gc())
parts$aoki <- to_universe(counts, "aoki")
meta$aoki <- data.table(cell = paste0("aoki_", colnames(counts)), dataset = "aoki", sample = cell_meta$StudyID,
                        chemistry = "3prime", label = cell_meta$final_2025_08_12, label_source = "published (Fig S2B)")


# ---- Step 3. Stewart: rebuild raw counts (B0) -------------------------------------------------------------
# The Stewart object holds log1p-normalised values. Undo the log (expm1), rescale each cell so it sums to its
# total_counts, round to whole counts, and drop the zeros.
stewart_obj <- readRDS(file.path(W, "input/stewart_filtered.rds"))
lognorm <- LayerData(stewart_obj, assay = "RNA", layer = "counts")
cell_meta <- stewart_obj@meta.data
rm(stewart_obj)
invisible(gc())
unlogged <- lognorm
unlogged@x <- expm1(unlogged@x)
raw <- unlogged %*% Diagonal(x = cell_meta$total_counts / colSums(unlogged))
raw@x <- round(raw@x)
raw <- drop0(raw)
stopifnot(all(abs(colSums(raw) - cell_meta$total_counts) <= 0.01 * cell_meta$total_counts))   # within 1% of total_counts
dimnames(raw) <- dimnames(lognorm)
parts$stewart <- to_universe(raw, "stewart")
meta$stewart <- data.table(cell = paste0("stewart_", colnames(lognorm)), dataset = "stewart", sample = cell_meta$channel,
                           chemistry = "3prime", label = cell_meta$final_2025_08_12, label_source = "published (Fig S2C)")


# ---- Step 4. 10x ---------------------------------------------------------------------------------------------
# The generic 10x "B_cell" label becomes "B_preGC" and is flagged ref_exclude (S2D convention only).
tenx_obj <- readRDS(file.path(W, "input/tenx_filtered.rds"))
counts <- LayerData(tenx_obj, assay = "RNA", layer = "counts")
cell_meta <- tenx_obj@meta.data
rm(tenx_obj)
invisible(gc())
parts$tenx <- to_universe(counts, "tenx")
tenx_label <- as.character(cell_meta$final_cluster2)
meta$tenx <- data.table(cell = paste0("tenx_", colnames(counts)), dataset = "tenx", sample = "10x_1", chemistry = "3prime",
                        label = ifelse(tenx_label == "B_cell", "B_preGC", tenx_label), label_source = "published (Fig S2A)",
                        ref_exclude = ifelse(tenx_label == "B_cell", "10x generic B_cell (S2D convention only)", NA_character_))


# ---- Step 5. Flex S1 (ex vivo Sample1) -------------------------------------------------------------------------
s1_obj <- readRDS(S1_OBJ)
s1_obj <- subset(s1_obj, subset = sample == "Sample1")
counts <- LayerData(s1_obj, assay = "RNA", layer = "counts")
cell_meta <- s1_obj@meta.data
rm(s1_obj)
invisible(gc())
stopifnot(ncol(counts) == 19435, all(counts@x == round(counts@x)))             # expected cell count; whole counts
parts$s1 <- to_universe(counts, "s1")

# S1 cell_type -> reference label. Types not in this map (Monocyte, Innate_T, Endothelial, IEG_stress) get label NA
# and a "provisional" label source (resolved at CP-A).
S1MAP <- c(CD4_T = "T_CD4", CD8_T = "T_CD8", Treg = "T_Treg", NK = "NK", Plasma = "Plasma", HRS = "HRS", Fibroblast = "CAF",
           cDC = "DC_Conventional", Macrophage = "Macrophage", B = "B_preGC")
s1_type <- as.character(cell_meta$cell_type)
meta$s1 <- data.table(cell = paste0("s1_", colnames(counts)), dataset = "s1", sample = "S1", chemistry = "Flex_v1.1",
                      label = unname(S1MAP[s1_type]),
                      label_source = ifelse(s1_type %in% names(S1MAP), "ex vivo f6 cell_type (mapped)", "provisional (resolve at CP-A)"),
                      s1_cell_type = s1_type, s1_cell_state = as.character(cell_meta$cell_state),
                      s1_state_matches_type = cell_meta$state_matches_type,
                      ref_exclude = ifelse(!cell_meta$state_matches_type, "S1 state_matches_type=FALSE", NA_character_))


# ---- Step 6. Merge into one Seurat object -----------------------------------------------------------------------
for (dataset in names(parts)) {
  colnames(parts[[dataset]]) <- meta[[dataset]]$cell
}
counts <- do.call(cbind, parts)                         # datasets side by side: lab, aoki, stewart, tenx, s1
rm(parts)
invisible(gc())
cell_table <- rbindlist(meta, fill = TRUE)              # columns missing in a dataset are filled with NA
stopifnot(identical(cell_table$cell, colnames(counts)))
rownames(counts) <- universe
obj <- CreateSeuratObject(counts, meta.data = data.frame(cell_table, row.names = cell_table$cell),
                          min.cells = 0, min.features = 0)     # no gene / cell filtering
log_("merged: %d genes x %d cells", nrow(obj), ncol(obj))
print(table(obj$dataset, obj$chemistry))


# ---- Step 7. Normalise and choose 3,000 HVGs ---------------------------------------------------------------------
obj <- NormalizeData(obj, verbose = FALSE)

## HVGs chosen within each dataset, ranked by how many datasets select them (as SelectIntegrationFeatures)
# 7a. 3,000 HVGs within each dataset (datasets in alphabetical order)
cells_by_dataset <- split(colnames(obj), obj$dataset)
hvg_by_dataset <- list()
for (dataset in names(cells_by_dataset)) {
  dataset_counts <- LayerData(obj, layer = "counts")[, cells_by_dataset[[dataset]]]
  dataset_obj <- FindVariableFeatures(CreateSeuratObject(dataset_counts), nfeatures = 3000, verbose = FALSE)
  hvg_by_dataset[[dataset]] <- VariableFeatures(dataset_obj)
}

# 7b. For every gene picked anywhere: in how many datasets it was picked, and its median rank in those datasets
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

# 7c. Top 3,000: most datasets first, ties broken by the better (lower) median rank
hvg <- candidate_genes[order(-n_datasets_per_gene, median_rank)][1:3000]
VariableFeatures(obj) <- hvg
fwrite(data.table(ensg = hvg, symbol = v47$symbol[match(hvg, v47$ensg)], n_datasets = as.integer(n_datasets_per_gene[hvg])),
       file.path(OUT, "hvg_3000.tsv"), sep = "\t")


# ---- Step 8. Scale, PCA, Harmony, clusters, UMAP -------------------------------------------------------------------
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
log_("integration done")


# ---- Step 9. Mixing: neighbour diversity (inverse Simpson over k = 30 harmony neighbours) ---------------------------
# For every cell: take its 30 nearest neighbours (Harmony space); the score is 1 / sum(p^2), where p = the share of
# the neighbours from each dataset (or chemistry). 1 = all neighbours from one group; higher = better mixed.
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


# ---- Step 10. Exports for CP-A (slim, no expression matrices) ----------------------------------------------------------
# 10a. per-cell table: metadata, all cluster columns, mixing scores, integrated and unintegrated UMAP coordinates
umap_harmony <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
cluster_cols <- grep("^RNA_snn_res", names(obj@meta.data), value = TRUE)
cells <- data.table(obj@meta.data[, c("dataset", "sample", "chemistry", "label", "label_source", "s1_cell_type", "s1_cell_state",
                                      "s1_state_matches_type", "ref_exclude", cluster_cols, "ilisi_dataset", "ilisi_chemistry")],
                    keep.rownames = "cell")
cells[, UMAP_1 := umap_harmony[, 1]]
cells[, UMAP_2 := umap_harmony[, 2]]
cells[, UMAPu_1 := umap_unintegrated[, 1]]
cells[, UMAPu_2 := umap_unintegrated[, 2]]
fwrite(cells, file.path(OUT, "cpA_cells.tsv.gz"), sep = "\t")

# 10b. dot-plot tables: for each marker, mean log-expression and % of cells expressing it, per group
markers <- c("PTPRC", "CD3E", "CD4", "CD8A", "FOXP3", "IL2RA", "MKI67", "TOP2A", "NKG7", "GNLY", "KLRD1", "TRDC", "TRGC1", "KIT",
             "IL7R", "MS4A1", "CD19", "TCL1A", "IGHD", "CD27", "BCL6", "AICDA", "RGS13", "MZB1", "JCHAIN", "SDC1", "CD68", "CD163",
             "C1QA", "LYZ", "CD14", "VCAN", "FCGR3A", "CLEC9A", "CD1C", "LILRA4", "IL3RA", "COL1A1", "COL1A2", "PDGFRA", "PECAM1",
             "VWF", "TNFRSF8", "IRF4", "CD274", "FOS", "JUN", "HSPA1A", "EGR1")
marker_genes <- v47[symbol %in% markers & ensg %in% rownames(obj)]          # in gene-table order
marker_data <- LayerData(obj, layer = "data")[marker_genes$ensg, ]          # log-normalised values

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
fwrite(dotplot_table(obj$RNA_snn_res.0.6), file.path(OUT, "cpA_dotplot_by_cluster_res0.6.tsv"), sep = "\t")
# by published label; S1 cells without a label are grouped by their provisional type ("S1:<type>")
label_groups <- ifelse(is.na(obj$label), paste0("S1:", obj$s1_cell_type), obj$label)
fwrite(dotplot_table(label_groups), file.path(OUT, "cpA_dotplot_by_label.tsv"), sep = "\t")

# 10c. save the object and the session info
saveRDS(obj, file.path(OUT, "cHL_merged_harmony_v2.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "cpA_sessionInfo.txt"))
log_("saved object and CP-A tables")
