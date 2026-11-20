# =============================================================================
# v2_C0_qc_scores.R — QC, lineage scores and contamination flags on round r2, before the sub-UMAP cleanup (Sherlock).
#
# What it adds to the round-2 object:
#  1. per-cell QC from each dataset's ORIGINAL counts (the merged object holds only the protein-coding Flex-panel universe):
#     nCount_orig, nFeature_orig, percent.mt, percent.ribo
#  2. lineage module scores on the merged object (T, NK, B, plasma, myeloid, DC, stromal, HRS) + IEG/stress score
#  3. major type per r2 cluster (author decision 2026-10-13)
#  4. contamination flag: a cell is flagged for lineage L when L is NOT one of its major type's own lineages and its L score is
#     at least the 25th percentile of L scores among cells published as L (e.g. HRS-labelled cells), i.e. as high as a genuine L cell.
#  5. protected cells (author 2026-10-13): cells published as T_gd or HRS (rare types) are never removed automatically;
#     their flags are reported separately.
# No doublet calling (author decision).
# Input : output/cHL_merged_harmony_v2_r2.rds, the five input objects (lab, Aoki, 10x, Stewart, ex vivo S1)
# Output: output/cHL_merged_v2_r2_qc.rds, output/cpB_cell_qc.tsv.gz, output/cpB_lineage_sets.tsv
# Run   : Rscript v2_C0_qc_scores.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})
W <- "${PROJECT_ROOT}/cHL_v2"             # project folder on Sherlock
OUT <- file.path(W, "output")
IN <- file.path(W, "input")
S1_OBJ <- "${PROJECT_ROOT}/scRNA/2026-10-04_cHL_ex_vivo/output/2026-10-04_chl_exvivo_final_seurat_obj_f6.rds"
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))

# log_(): print a time-stamped message
log_ <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

obj <- readRDS(file.path(OUT, "cHL_merged_harmony_v2_r2.rds"))
cells <- colnames(obj)


# ---- Step 1. QC from the original counts ----------------------------------------------------------------------
# qc_from_counts(): per cell, total UMIs, genes detected, % mitochondrial (MT-) and % ribosomal (RPS/RPL) UMIs.
#   prefix = the dataset prefix the merged object puts in front of the barcode (e.g. "lab_").
qc_from_counts <- function(m, prefix) {
  is_mt   <- grepl("^MT-", rownames(m))
  is_ribo <- grepl("^RP[SL][0-9]", rownames(m))
  total   <- Matrix::colSums(m)
  data.table(cell = paste0(prefix, colnames(m)),
             nCount_orig = total,
             nFeature_orig = Matrix::colSums(m > 0),
             percent.mt = 100 * Matrix::colSums(m[is_mt, , drop = FALSE]) / pmax(total, 1),     # pmax: avoid / 0
             percent.ribo = 100 * Matrix::colSums(m[is_ribo, , drop = FALSE]) / pmax(total, 1))
}

qc <- list()
# (gc() is left visible as in the original script, so its memory table is printed to the log)
input_obj <- readRDS(file.path(IN, "lab_scRNA_seurat.rds"))
qc$lab <- qc_from_counts(LayerData(input_obj, assay = "RNA", layer = "counts"), "lab_")
rm(input_obj)
gc()
input_obj <- readRDS(file.path(IN, "aoki_filtered.rds"))
qc$aoki <- qc_from_counts(LayerData(input_obj, assay = "originalexp", layer = "counts"), "aoki_")
rm(input_obj)
gc()
input_obj <- readRDS(file.path(IN, "tenx_filtered.rds"))
qc$tenx <- qc_from_counts(LayerData(input_obj, assay = "RNA", layer = "counts"), "tenx_")
rm(input_obj)
gc()

# Stewart has ENSG rows and log1p values: rebuild raw counts as in step B0, then name the rows by v47 symbol
# so that MT- / RPS / RPL genes can be found
input_obj <- readRDS(file.path(IN, "stewart_filtered.rds"))
lognorm <- LayerData(input_obj, assay = "RNA", layer = "counts")
cell_meta <- input_obj@meta.data
rm(input_obj)
gc()
unlogged <- lognorm
unlogged@x <- expm1(unlogged@x)
raw <- unlogged %*% Diagonal(x = cell_meta$total_counts / colSums(unlogged))
raw@x <- round(raw@x)
raw <- drop0(raw)
dimnames(raw) <- dimnames(lognorm)
rownames(raw) <- v47$symbol[match(sub("\\.[0-9]+$", "", rownames(raw)), v47$ensg)]
rownames(raw)[is.na(rownames(raw))] <- "NA_gene"                    # ENSG without a v47 symbol
qc$stewart <- qc_from_counts(raw, "stewart_")
rm(raw, unlogged, lognorm)
gc()

# ex vivo S1: Sample1 only
input_obj <- readRDS(S1_OBJ)
input_obj <- subset(input_obj, subset = sample == "Sample1")
qc$s1 <- qc_from_counts(LayerData(input_obj, assay = "RNA", layer = "counts"), "s1_")
rm(input_obj)
gc()

# put the QC rows in the object's cell order and attach them
qc <- rbindlist(qc)
qc <- qc[match(cells, cell)]
stopifnot(!anyNA(qc$cell))                                          # every cell found
for (col in c("nCount_orig", "nFeature_orig", "percent.mt", "percent.ribo")) {
  obj[[col]] <- qc[[col]]
}
log_("QC attached for %d cells (MT genes: Flex panels carry no MT- probes, so percent.mt is ~0 for Flex cells)", length(cells))


# ---- Step 2. Lineage + stress module scores -------------------------------------------------------------------------
# sym2ensg(): gene symbols -> ENSG IDs that are in the object (symbols not found are dropped)
sym2ensg <- function(symbols) {
  ensg <- v47$ensg[match(symbols, v47$symbol)]
  ensg[!is.na(ensg) & ensg %in% rownames(obj)]
}
LIN <- list(
  T       = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD6", "LCK", "ZAP70", "CD247"),
  NK      = c("NKG7", "GNLY", "KLRD1", "PRF1", "GZMB", "KLRF1", "NCAM1", "FCGR3A"),
  B       = c("MS4A1", "CD19", "CD79A", "CD79B", "BANK1", "PAX5", "CD22", "FCRL1"),
  Plasma  = c("MZB1", "JCHAIN", "XBP1", "PRDM1", "SDC1", "TNFRSF17", "DERL3"),
  Myeloid = c("LYZ", "CD68", "CD14", "CD163", "C1QA", "C1QB", "FCER1G", "TYROBP", "AIF1", "CSF1R"),
  DC      = c("CLEC9A", "XCR1", "CD1C", "FCER1A", "CLEC10A", "LILRA4", "IL3RA", "LAMP3"),
  Stromal = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "PDGFRB", "FAP"),
  HRS     = c("TNFRSF8", "CCL17", "CCL22", "IL13", "IRF4", "BATF3", "FSCN1", "CD274"),
  Stress  = c("FOS", "FOSB", "JUN", "JUNB", "EGR1", "IER2", "HSPA1A", "HSPA1B", "DNAJB1", "ATF3"))
sets <- list()
for (set_name in names(LIN)) {
  sets[[set_name]] <- sym2ensg(LIN[[set_name]])
}
print(lengths(sets))                                                # genes per set found in the universe

# AddModuleScore names the scores LINSCORE1, LINSCORE2, ... (seed 42); rename them to <set>_score (e.g. T_score)
obj <- AddModuleScore(obj, features = sets, name = "LINSCORE", seed = 42)
for (i in seq_along(sets)) {
  new_name <- paste0(names(sets)[i], "_score")
  old_name <- paste0("LINSCORE", i)
  obj[[new_name]] <- obj[[old_name]]
  obj[[old_name]] <- NULL
}


# ---- Step 3. Major type per r2 cluster (author decision 2026-10-13) ---------------------------------------------------
MAJOR <- c(`1` = "T_CD4_Treg", `3` = "T_CD4_Treg", `7` = "T_CD4_Treg",
           `4` = "T_CD8_NK_gd", `9` = "T_CD8_NK_gd", `15` = "T_CD8_NK_gd",
           `8` = "T_MKI67",
           `0` = "B_preGC_Memory", `2` = "B_preGC_Memory", `5` = "B_preGC_Memory",
           `11` = "Plasma_BGC_HRS", `13` = "Plasma_BGC_HRS", `18` = "Plasma_BGC_HRS", `16` = "Plasma_BGC_HRS",
           `6` = "Myeloid", `10` = "Myeloid", `12` = "Myeloid", `17` = "Myeloid",
           `14` = "Stromal")
obj$r2_cluster <- as.character(obj$RNA_snn_res.0.6)
obj$major_type <- unname(MAJOR[obj$r2_cluster])
stopifnot(!anyNA(obj$major_type))                                   # every r2 cluster has a major type
print(table(obj$major_type, obj$dataset))


# ---- Step 4. Contamination flags --------------------------------------------------------------------------------------
# the lineages each major type is allowed to express ("own" lineages)
OWN <- list(T_CD4_Treg = "T", T_CD8_NK_gd = c("T", "NK"), T_MKI67 = c("T", "NK"), B_preGC_Memory = "B",
            Plasma_BGC_HRS = c("B", "Plasma", "HRS"), Myeloid = c("Myeloid", "DC"), Stromal = "Stromal")
## threshold = 25th percentile of L among the cells whose PUBLISHED LABEL is L (not the whole owning major type: e.g. the
## HRS owner Plasma_BGC_HRS is mostly plasma/B-GC, so a major-type quantile sits at background and flagged ~all cells)
OWNER_LABELS <- list(T = c("T_CD4", "T_CD8", "T_gd", "T_MKI67", "T_Treg"), NK = "NK", B = c("B_preGC", "B_Memory", "B_GC"),
                     Plasma = "Plasma", Myeloid = "Macrophage", DC = c("DC_Conventional", "DC_Plasmacytoid"), Stromal = "CAF", HRS = "HRS")
lineages <- setdiff(names(LIN), "Stress")
stopifnot(setequal(lineages, names(OWNER_LABELS)))
md <- data.table(obj@meta.data, keep.rownames = "cell")             # per-cell table

# 4a. threshold per lineage
thr <- numeric(0)
for (lineage in lineages) {
  owner_cells <- md[label %in% OWNER_LABELS[[lineage]]]
  owner_scores <- owner_cells[[paste0(lineage, "_score")]]
  thr[lineage] <- unname(quantile(owner_scores, 0.25))   # quantile() names its value "25%": unname so thr["T"] works
}
stopifnot(identical(names(thr), lineages), !anyNA(thr))
print(round(thr, 3))

# 4b. flag matrix (cells x lineages): TRUE when the lineage is foreign to the cell's major type AND the score >= threshold
flag <- matrix(FALSE, nrow = nrow(md), ncol = length(lineages), dimnames = list(NULL, lineages))
for (lineage in lineages) {
  owning_major_types <- c()                               # major types that list this lineage as their own
  for (major in names(OWN)) {
    if (lineage %in% OWN[[major]]) owning_major_types <- c(owning_major_types, major)
  }
  foreign <- !(md$major_type %in% owning_major_types)
  flag[, lineage] <- foreign & md[[paste0(lineage, "_score")]] >= thr[lineage]
}
stopifnot(!anyNA(flag))

# 4c. per cell: the flagged lineages as text (e.g. "B,Myeloid"; "" = none) and contaminated = any flag
flagged_text <- character(nrow(flag))
for (i in seq_len(nrow(flag))) {
  flagged_text[i] <- paste(lineages[flag[i, ]], collapse = ",")
}
md[, contam_lineages := flagged_text]
md[, contaminated := contam_lineages != ""]

## author 2026-10-13: T_gd and HRS are rare -> protected: never removed automatically; their flags are reported separately
md[, protected := label %in% c("T_gd", "HRS")]
obj$protected <- md$protected
cat("protected cells (T_gd + HRS):", sum(md$protected), "; of those flagged:", sum(md$protected & md$contaminated), "\n")
obj$contam_lineages <- md$contam_lineages
obj$contaminated <- md$contaminated

# 4d. print: % flagged per major type, and the 15 most frequent (major type, flagged lineages) combinations
print(md[, .(cells = .N, pct_flagged = round(100 * mean(contaminated), 1)), by = major_type])
flag_combinations <- md[contaminated == TRUE, .N, by = .(major_type, contam_lineages)]
flag_combinations <- flag_combinations[order(-N)]
print(head(flag_combinations, 15))


# ---- Step 5. Save ------------------------------------------------------------------------------------------------------
qc_cols <- c("cell", "dataset", "sample", "label", "r2_cluster", "major_type", "nCount_orig", "nFeature_orig", "percent.mt",
             "percent.ribo", paste0(names(LIN), "_score"), "contaminated", "contam_lineages", "protected")
fwrite(md[, qc_cols, with = FALSE], file.path(OUT, "cpB_cell_qc.tsv.gz"), sep = "\t")

# lineage gene sets: genes (comma-separated), how many are in the universe, and the flag threshold (none for Stress)
set_genes <- character(0)
for (set_name in names(LIN)) {
  set_genes[set_name] <- paste(LIN[[set_name]], collapse = ",")
}
fwrite(data.table(lineage = names(LIN), genes = set_genes, n_in_universe = lengths(sets),
                  flag_threshold = c(round(thr, 4), NA)),
       file.path(OUT, "cpB_lineage_sets.tsv"), sep = "\t")
saveRDS(obj, file.path(OUT, "cHL_merged_v2_r2_qc.rds"))
log_("saved r2 object with QC, lineage scores, major types and contamination flags")
