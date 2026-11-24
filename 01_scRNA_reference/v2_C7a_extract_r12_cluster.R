# =============================================================================
# v2_C7a_extract_r12_cluster.R
#
# PURPOSE : Cut one round-12 cluster out of the 1.4 GB object (counts + metadata only) so the sub-UMAP (v2_C7) can run on a
#           laptop (author 2026-10-20).
# INPUT   : output/cHL_merged_harmony_v2_r12.rds (on OAK); env CLUSTER = res-0.6 cluster to keep (default 8)
# OUTPUT  : output/subumap/r12_cl<CLUSTER>_input.rds
# RUN     : on Sherlock, CLUSTER=8 Rscript v2_C7a_extract_r12_cluster.R
# =============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})
W  <- "${PROJECT_ROOT}/cHL_v2"      # project root on OAK
CL <- Sys.getenv("CLUSTER", "8")

full_obj <- readRDS(file.path(W, "output/cHL_merged_harmony_v2_r12.rds"))

# Step 1. Cells of the chosen round-12 cluster (res 0.6)
keep <- colnames(full_obj)[as.character(full_obj$RNA_snn_res.0.6) == CL]

# Step 2. New small object: raw counts + all metadata of those cells (no reductions)
cluster_obj <- CreateSeuratObject(LayerData(full_obj, assay = "RNA", layer = "counts")[, keep],
                                  meta.data = full_obj@meta.data[keep, ])
saveRDS(cluster_obj, file.path(W, "output/subumap", sprintf("r12_cl%s_input.rds", CL)))
cat("saved", ncol(cluster_obj), "cells\n")
