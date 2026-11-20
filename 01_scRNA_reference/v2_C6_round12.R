# =============================================================================
# v2_C6_round12.R — round 12 of the merged cHL object (runs on Sherlock).
#
# PURPOSE (author 2026-10-16)
#   Take round 11 (r11_rolling), apply the results of the cluster 4/8/16 (CD8 / NK / ILC / T_gd) sub-UMAP and re-integrate
#   with Harmony theta 4/4/2:
#   1. Remove the 945 cells dropped in the sub-UMAP cleanup rounds (r2 253 heat-shock CD8; r3 104 MT-dominated; r4 588
#      generic weak + CD4 T cells in the CD8 cluster); file output/stepC/inputs/2026-10-16_cl4_8_16_removed_cells_all_rounds.tsv.gz.
#   2. Final labels = cluster annotation of r11 (cluster type = majority published label of the cluster's non-S1 cells;
#      protected HRS / T_gd / ILC keep their own published label), replaced for the 13,619 sub-UMAP cells by the sub-UMAP
#      labels (subcluster 14 -> T_gd; subcluster 3 sub-subclusters 2/4/6 -> T_gd); file .../2026-10-16_cl4_8_16_final_labels.tsv.gz.
#   3. Recluster the whole object with the rounds 1-11 recipe (per-dataset HVG 3,000 -> PCA 50 -> Harmony -> Louvain res
#      0.3 / 0.6 / 1.0 -> UMAP), but Harmony on dataset + chemistry + sample with theta 4 / 4 / 2 (rounds 1-11: dataset +
#      chemistry, default theta 2; 4/4/2 chosen by the author from the sub-UMAP theta test, v2_C5g).
#
# INPUT  : output/cHL_merged_harmony_v2_r11_rolling.rds, the two sub-UMAP files above, ref/gene_table_v47.tsv
# OUTPUT : output/cHL_merged_harmony_v2_r12.rds, output/cpA12_* (same layout as cpA11rolling_, + final_label)
# RUN    : on Sherlock (sbatch v2_run.sh v2_C6_round12.R)
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
IN  <- file.path(OUT, "stepC/inputs")
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))

# Print a time-stamped log line
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

THETA     <- c(dataset = 4, chemistry = 4, sample = 2)     # Harmony theta per batch variable (author, from v2_C5g)
PROTECTED <- c("HRS", "T_gd", "ILC")                        # rare types that keep their own published label
SEED      <- 112L                                           # seed for PCA, Harmony and clustering in this round


## Step 1. Round-11 cluster annotation ------------------------------------------------------------------------------------
obj <- readRDS(file.path(OUT, "cHL_merged_harmony_v2_r11_rolling.rds"))
all_meta   <- data.table(obj@meta.data, keep.rownames = "cell")
cell_table <- all_meta[, .(cell, dataset, label, r11 = as.character(RNA_snn_res.0.6))]
rm(all_meta)

# Cluster type = the most common published label among the cluster's non-S1 labelled cells
label_counts <- cell_table[dataset != "s1" & !is.na(label), .N, by = .(r11, label)]
label_counts <- label_counts[order(r11, -N)]                       # most common label first within each cluster
ctype <- label_counts[, .(cluster_type = label[1]), by = r11]
cell_table <- merge(cell_table, ctype, by = "r11")
cell_table <- cell_table[match(colnames(obj), cell)]                # back to object order

# Final label = cluster type; protected rare types (non-S1) keep their own published label
cell_table[, final_label := fifelse(dataset != "s1" & label %in% PROTECTED, label, cluster_type)]
cell_table[, label_source := "r11 cluster annotation (protected types: own label)"]


## Step 2. Apply the cluster 4/8/16 sub-UMAP: its labels, and its removed cells ----------------------------------------------
subumap_removed <- fread(cmd = sprintf("gzip -dc %s", file.path(IN, "2026-10-16_cl4_8_16_removed_cells_all_rounds.tsv.gz")))
subumap_labels  <- fread(cmd = sprintf("gzip -dc %s", file.path(IN, "2026-10-16_cl4_8_16_final_labels.tsv.gz")))
stopifnot(nrow(subumap_removed) == 945,
          nrow(subumap_labels) == 13619,
          all(subumap_removed$cell %in% cell_table$cell),
          all(subumap_labels$cell %in% cell_table$cell),
          !any(subumap_labels$cell %in% subumap_removed$cell))
# Join on cell: the sub-UMAP label (and its source) replaces the cluster label
cell_table[subumap_labels, on = "cell", `:=`(final_label  = i.final_label,
                                             label_source = paste0("cluster 4/8/16 sub-UMAP: ", i.source))]
obj$final_label        <- cell_table$final_label
obj$final_label_source <- cell_table$label_source
obj$r11_cluster        <- cell_table$r11
fwrite(subumap_removed, file.path(OUT, "cpA12_removed_cells.tsv.gz"), sep = "\t")
obj <- subset(obj, cells = setdiff(colnames(obj), subumap_removed$cell))
log_("r11: %d cells; removed %d; r12 input %d cells", nrow(cell_table), nrow(subumap_removed), ncol(obj))
print(table(final_label = obj$final_label, dataset = obj$dataset))

# Keep the r11 clusters under a new name; drop the old clustering columns, reductions, graphs and neighbours
obj$r11_cluster_res0.6 <- obj$r11_cluster
old_cluster_columns <- grep("^RNA_snn_res|^seurat_clusters$", names(obj@meta.data), value = TRUE)
for (column in old_cluster_columns) {
  obj[[column]] <- NULL
}
obj@reductions <- list()
obj@graphs     <- list()
obj@neighbors  <- list()


## Step 3. Recluster: rounds 1-11 recipe, Harmony dataset + chemistry + sample, theta 4/4/2 ----------------------------------
# 3a. 3,000 variable genes per dataset
cells_by_dataset <- split(colnames(obj), obj$dataset)
hv <- list()                                                     # one HVG vector per dataset
for (ds in names(cells_by_dataset)) {
  ds_obj <- CreateSeuratObject(LayerData(obj, assay = "RNA", layer = "counts")[, cells_by_dataset[[ds]]])
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
DefaultAssay(obj) <- "RNA"
VariableFeatures(obj) <- hvg
fwrite(data.table(ensg = hvg, symbol = v47$symbol[match(hvg, v47$ensg)], n_datasets = as.integer(rk[hvg])),
       file.path(OUT, "cpA12_hvg_3000.tsv"), sep = "\t")

# 3c. Scale, PCA, unintegrated UMAP, Harmony, Louvain clustering, integrated UMAP
obj <- ScaleData(obj, features = hvg, verbose = FALSE)
obj <- RunPCA(obj, features = hvg, npcs = 50, seed.use = SEED, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 42, reduction.name = "umap_unintegrated", verbose = FALSE)
set.seed(SEED)                                                   # Harmony uses R's random numbers
obj <- RunHarmony(obj, group.by.vars = names(THETA), theta = unname(THETA), reduction.use = "pca", dims.use = 1:50, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
for (r in c(0.3, 0.6, 1.0)) {
  obj <- FindClusters(obj, resolution = r, random.seed = SEED, verbose = FALSE)
}
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, seed.use = 42, verbose = FALSE)

# 3d. Mixing score (iLISI): per cell, inverse Simpson index of dataset / chemistry among its 30 nearest neighbours
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
log_("r12: %d cells, %d clusters at res 0.6; median iLISI dataset %.2f (S1 cells %.2f)", ncol(obj), length(unique(obj$RNA_snn_res.0.6)),
     median(obj$ilisi_dataset), median(obj$ilisi_dataset[obj$dataset == "s1"]))


## Step 4. Exports (cpA11rolling_ layout + final_label) -------------------------------------------------------------------
umap_integrated   <- Embeddings(obj, "umap")
umap_unintegrated <- Embeddings(obj, "umap_unintegrated")
cl <- grep("^RNA_snn_res", names(obj@meta.data), value = TRUE)
keep_cols <- intersect(c("dataset", "sample", "chemistry", "label", "label_source", "s1_cell_type", "s1_cell_state", "s1_state_matches_type", "ref_exclude",
                         "r1_cluster_res0.6", "r2_cluster_res0.6", "major_type", "r2_subcluster", "r3_cluster", "r11_cluster_res0.6", "cell_type", "lineage",
                         "final_label", "final_label_source", cl, "ilisi_dataset", "ilisi_chemistry", "nCount_RNA", "nFeature_RNA", "nCount_orig",
                         "nFeature_orig", "percent.mt", "contaminated", "contam_lineages", "protected"), names(obj@meta.data))
cells <- data.table(obj@meta.data[, keep_cols], keep.rownames = "cell")
cells[, `:=`(UMAP_1  = umap_integrated[, 1],   UMAP_2  = umap_integrated[, 2],
             UMAPu_1 = umap_unintegrated[, 1], UMAPu_2 = umap_unintegrated[, 2])]
fwrite(cells, file.path(OUT, "cpA12_cells.tsv.gz"), sep = "\t")

# Each res-0.6 cluster: cells per final label and their % of the cluster
chk <- cells[, .N, by = .(cluster = RNA_snn_res.0.6, final_label)]
chk[, pct := round(100 * N / sum(N), 1), by = cluster]
chk <- chk[order(cluster, -N)]
fwrite(chk, file.path(OUT, "cpA12_cluster_by_final_label.tsv"), sep = "\t")
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
fwrite(dot(obj$RNA_snn_res.0.6), file.path(OUT, "cpA12_dotplot_by_cluster_res0.6.tsv"), sep = "\t")
fwrite(dot(obj$final_label), file.path(OUT, "cpA12_dotplot_by_final_label.tsv"), sep = "\t")
# by published label; S1 cells (no published label) by "S1:<ex vivo type>"
fwrite(dot(ifelse(is.na(obj$label), paste0("S1:", obj$s1_cell_type), obj$label)), file.path(OUT, "cpA12_dotplot_by_label.tsv"), sep = "\t")
saveRDS(obj, file.path(OUT, "cHL_merged_harmony_v2_r12.rds"))
log_("saved r12 object and cpA12_ tables")
