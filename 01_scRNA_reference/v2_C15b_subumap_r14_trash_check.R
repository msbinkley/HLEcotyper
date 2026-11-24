# =============================================================================
# v2_C15b_subumap_r14_trash_check.R — trash-cell check of the res-1.2 subclusters of the round-14 T_CD8 / NK / ILC
# sub-UMAP (v2_C15; author 2026-10-25: "are there trash cells in subclusters, res 1.2").
#
# Per subcluster (res 1.2), using the sub-UMAP object (log-normalised RNA):
#   - QC: median UMI, genes, mito %
#   - lineage scores = mean log-expression of a small marker panel per foreign lineage (B, plasma, myeloid, HRS,
#     CAF / stroma, erythroid, platelet) and of the own lineages (T, NK, ILC)
#   - % of cells expressing the key foreign markers (a doublet / contamination sign is a high %)
#   - heat-shock / dissociation-stress score (as the v2 lineage cleanup, v2_C1)
# Nothing is removed. The author decides from the table and the plot.
#
# Run on Sherlock:  sbatch v2_run.sh v2_C15b_subumap_r14_trash_check.R
# Output: output/subumap/2026-10-25_r14_sub_cd8_nk_ilc_tgd_trash_check.tsv and _trash_check.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

# Step 1. The sub-UMAP object and the output names
work_dir <- "${PROJECT_ROOT}/cHL_v2"
out_dir  <- file.path(work_dir, "output/subumap")
out_file <- function(name) file.path(out_dir, paste0("2026-10-25_r14_sub_cd8_nk_ilc_tgd_", name))
sub <- readRDS(out_file("obj.rds"))
subcluster <- as.character(sub$RNA_snn_res.1.2)

# Step 2. Marker panels (only genes present in the object are used)
panels <- list(
  T_cell     = c("CD3D", "CD3E", "CD3G", "CD2", "TRAC"),
  NK         = c("KLRF1", "NCAM1", "FCGR3A", "NCR1", "KLRD1"),
  ILC        = c("KIT", "IL7R", "RORC", "IL1R1", "IL23R"),
  B_cell     = c("CD79A", "CD79B", "MS4A1", "CD19", "BANK1"),
  Plasma     = c("JCHAIN", "MZB1", "IGKC", "IGHG1", "XBP1"),
  Myeloid    = c("LYZ", "CD14", "CD68", "C1QA", "AIF1", "CST3", "S100A8"),
  HRS        = c("TNFRSF8", "IRF4", "CCL17", "CD274", "FSCN1"),
  CAF        = c("COL1A1", "COL3A1", "DCN", "LUM", "PDGFRB"),
  Erythroid  = c("HBB", "HBA1", "HBA2"),
  Platelet   = c("PPBP", "PF4"),
  HeatShock  = c("HSPA1A", "HSPA1B", "HSPA6", "HSPH1", "DNAJB1", "HSP90AA1"),
  Dissociation = c("FOS", "FOSB", "JUN", "ATF3", "EGR1", "IER2", "DUSP1", "NR4A1"),
  Metallothionein = c("MT1E", "MT1F", "MT1X", "MT2A")
)
key_foreign <- c("CD79A", "MS4A1", "JCHAIN", "LYZ", "CD14", "C1QA", "AIF1", "COL1A1", "HBB", "PPBP", "TNFRSF8")

# Step 3. Per subcluster: QC medians, panel scores (mean log-expression over the panel's genes, averaged over cells)
#         and % of cells expressing each key foreign marker
expr <- LayerData(sub, assay = "RNA", layer = "data")
rows <- list()
for (this_subcluster in sort(unique(as.integer(subcluster)))) {
  in_subcluster <- subcluster == as.character(this_subcluster)
  row <- data.table(
    subcluster  = this_subcluster,
    cells       = sum(in_subcluster),
    T_gd        = sum(sub$final_label[in_subcluster] == "T_gd"),
    median_umi  = median(sub$nCount_RNA[in_subcluster]),
    median_genes = median(sub$nFeature_RNA[in_subcluster]))
  for (panel_name in names(panels)) {
    genes <- intersect(panels[[panel_name]], rownames(expr))
    row[[paste0("score_", panel_name)]] <- round(mean(Matrix::colMeans(expr[genes, in_subcluster, drop = FALSE])), 3)
  }
  for (gene in intersect(key_foreign, rownames(expr))) {
    row[[paste0("pct_", gene)]] <- round(100 * mean(expr[gene, in_subcluster] > 0), 1)
  }
  rows[[length(rows) + 1]] <- row
}
trash_table <- rbindlist(rows)
fwrite(trash_table, out_file("trash_check.tsv"), sep = "\t")
print(trash_table)

# Step 4. Plot: panel scores per subcluster (dot plot of the marker panels) and QC violins
panel_genes <- unique(unlist(panels))
panel_genes <- panel_genes[panel_genes %in% rownames(sub)]
pdf(out_file("trash_check.pdf"), width = 18, height = 7)
print(DotPlot(sub, features = panel_genes, group.by = "RNA_snn_res.1.2") + RotatedAxis() +
      ggtitle("Lineage / stress panels by subcluster (res 1.2)"))
print(VlnPlot(sub, features = c("nCount_RNA", "nFeature_RNA"), group.by = "RNA_snn_res.1.2", pt.size = 0, log = TRUE))
dev.off()
cat("done\n")
