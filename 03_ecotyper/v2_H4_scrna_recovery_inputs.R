# =============================================================================
# v2_H4_scrna_recovery_inputs.R — merged scRNA input of ONE round (16 or 17) for the EcoTyper scRNA recovery of the
# hit models (Fig 3e); author 2026-10-28 (per-hit figure sets; Fig 3e with EcoTyper -s 2500 downsampling).
#
# Same as v2_G4_scrna_recovery_inputs.R (round 14), with the round as a parameter:
#   count      = raw counts of every cell of the round (v47 symbols), Cycling cells left out
#   annotation = ID, CellType (the 13-type model names: B_postGC = B_GC + B_Memory, NK = NK + ILC), Sample, Dataset
# Round 17 also gets a second annotation for the DC-pooled model (r17 DC hit): DC_Conventional + DC_Plasmacytoid = DC
# (as v2_D17_pool12_dc.R pools the bulk fractions). The count file is shared by both annotations.
#
# Run on Sherlock (in the job of v2_H4_round_extras.sh):  ROUND=r16 Rscript v2_H4_scrna_recovery_inputs.R
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
})

# Step 1. Round, files, expected cell numbers
ROUND <- Sys.getenv("ROUND")
stopifnot(ROUND %in% c("r16", "r17"))
V2_DIR      <- "${PROJECT_ROOT}/cHL_v2"
OBJECT_FILE <- file.path(V2_DIR, "output", sprintf("cHL_merged_harmony_v2_%s.rds", ROUND))
CELLS_FILE  <- file.path(V2_DIR, "output", sprintf("2026-10-26_cpA%s_cells.tsv.gz", sub("r", "", ROUND)))   # cell, dataset,
                                                                                                # sample, final_label
GENE_TABLE  <- file.path(V2_DIR, "ref/gene_table_v47.tsv")                                      # ensg -> v47 symbol
OUTPUT_DIR  <- "${PROJECT_ROOT}/Ecotyper/datasets"
COUNT_FILE  <- file.path(OUTPUT_DIR, sprintf("2026-10-28_cHL_%s_scRNA_merged_count.txt", ROUND))
ANNOT_FILE  <- file.path(OUTPUT_DIR, sprintf("2026-10-28_cHL_%s_scRNA_merged_annotation.txt", ROUND))
ANNOT_DC_FILE <- file.path(OUTPUT_DIR, sprintf("2026-10-28_cHL_%s_scRNA_merged_annotation_DC.txt", ROUND))
EXPECTED_CELLS <- c(r16 = 110043, r17 = 109952)      # v2_C21 (round 16) and v2_C23 (round 17)
GENE_BLOCK  <- 1000                                  # genes written per block (memory cap)

# Step 2. The round's cells (dataset, sample, final label) and their raw counts (LayerData, as v2_C4)
cell_table <- fread(cmd = sprintf("gzip -dc %s", CELLS_FILE), select = c("cell", "dataset", "sample", "final_label"))
stopifnot(nrow(cell_table) == EXPECTED_CELLS[[ROUND]])
object <- readRDS(OBJECT_FILE)
counts <- LayerData(object, assay = "RNA", layer = "counts")
stopifnot(all(cell_table$cell %in% colnames(counts)))
counts <- counts[, cell_table$cell]
rm(object)
invisible(gc())

# Step 3. Rows: Ensembl ID -> v47 symbol (every gene must map, and symbols must be unique)
gene_table <- fread(GENE_TABLE)
symbols <- gene_table$symbol[match(rownames(counts), gene_table$ensg)]
stopifnot(!anyNA(symbols), !anyDuplicated(symbols))
rownames(counts) <- symbols

# Step 4. Labels of the 13-type models (as v2_G4); Cycling cells left out
cell_types <- cell_table$final_label
cell_types[cell_types %in% c("B_GC", "B_Memory")] <- "B_postGC"
cell_types[cell_types == "ILC"] <- "NK"
keep <- cell_types != "Cycling"
cat(sprintf("%s: cells kept %d of %d (Cycling left out)\n", ROUND, sum(keep), length(keep)))
print(table(cell_types[keep], cell_table$dataset[keep]))

# Step 5. Annotation (cell IDs with "-" -> ".", as EcoTyper reads them)
cell_ids <- gsub("-", ".", cell_table$cell[keep])
annotation <- data.table(ID = cell_ids, CellType = cell_types[keep], Sample = cell_table$sample[keep],
                         Dataset = cell_table$dataset[keep])
fwrite(annotation, ANNOT_FILE, sep = "\t")

# Step 6. Round 17 only: the DC-pooled annotation (cDC + pDC = DC) for the r17 DC model
if (ROUND == "r17") {
  annotation_dc <- copy(annotation)
  annotation_dc[CellType %in% c("DC_Conventional", "DC_Plasmacytoid"), CellType := "DC"]
  fwrite(annotation_dc, ANNOT_DC_FILE, sep = "\t")
  cat("DC-pooled annotation:", sum(annotation_dc$CellType == "DC"), "DC cells\n")
}

# Step 7. Counts: header line, then blocks of GENE_BLOCK genes (dense only within a block)
counts <- counts[, keep]
writeLines(paste(c("Gene", cell_ids), collapse = "\t"), COUNT_FILE)
for (block_start in seq(1, nrow(counts), by = GENE_BLOCK)) {
  block_rows <- block_start:min(block_start + GENE_BLOCK - 1, nrow(counts))
  block <- as.matrix(counts[block_rows, , drop = FALSE])
  fwrite(data.table(Gene = rownames(block), block), COUNT_FILE, sep = "\t", append = TRUE, col.names = FALSE)
}
cat(sprintf("wrote %s (%d genes x %d cells) and %s\n", COUNT_FILE, nrow(counts), ncol(counts), ANNOT_FILE))
