# =============================================================================
# v3_I1_cytospace_inputs_r16.R — CytoSPACE scRNA inputs from ROUND 16 (author 2026-10-30: "do both" [re-run cfRNA
# and CytoSPACE with the round-16 reference]; "use all 110,043 r16 cells").
#
# v1 (Cytospace/2025-07-12_cHL_cytospace.sh) gave CytoSPACE all genes x all cells of the v1 reference (44,769 cells,
# raw counts, 16 labels). Here: the same layout from round 16:
#   - counts: the raw UMI counts on the gene universe saved by the round-16 reference run (v2_C4 -> output/r16/stepC/
#     cHL_v2_r10_counts_universe.rds; the same cells and counts that built the round-16 CIBERSORTx reference)
#   - labels: the round-16 final labels (16 types; "Cycling" is v1's "T_MKI67+")
# Outputs (OAK cHL_v2/output/r16/stepI/cytospace_input/):
#   2026-10-30_r16_scRNA_data.txt         GENES + one column per cell, raw counts (tab, as v1)
#   2026-10-30_r16_cell_type_labels.txt   "Cell IDs" <tab> "CellType" (as v1)
# Run on Sherlock:  sbatch v3_I1_cytospace_inputs_r16.sh
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

# Step 1. Files
V2_DIR       <- "${PROJECT_ROOT}/cHL_v2"
COUNTS_FILE  <- file.path(V2_DIR, "output/r16/stepC/cHL_v2_r10_counts_universe.rds")
LABELS_FILE  <- file.path(V2_DIR, "output/stepC/inputs/2026-10-26_round16_labels.tsv.gz")
OUTPUT_DIR   <- file.path(V2_DIR, "output/r16/stepI/cytospace_input")
DATA_FILE    <- file.path(OUTPUT_DIR, "2026-10-30_r16_scRNA_data.txt")
LABEL_FILE   <- file.path(OUTPUT_DIR, "2026-10-30_r16_cell_type_labels.txt")
GENES_PER_CHUNK <- 500                       # genes written per block (500 genes x 110,043 cells fits easily in memory)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Step 2. Counts and labels; keep the cells in the order of the counts matrix
saved  <- readRDS(COUNTS_FILE)
counts <- saved$counts
rm(saved)
invisible(gc())
labels <- fread(cmd = sprintf("gzip -dc '%s'", LABELS_FILE), select = c("cell", "final_label"))
cat(sprintf("counts: %d genes x %d cells; labels: %d cells\n", nrow(counts), ncol(counts), nrow(labels)))
stopifnot(setequal(colnames(counts), labels$cell))         # the reference cells are exactly the round-16 cells
labels <- labels[match(colnames(counts), labels$cell)]
stopifnot(nrow(labels) == 110043, !anyNA(labels$final_label))
print(labels[, .N, by = final_label][order(-N)])

# Step 3. Labels file (header as v1)
label_table <- data.table(`Cell IDs` = labels$cell, CellType = labels$final_label)
fwrite(label_table, LABEL_FILE, sep = "\t")

# Step 4. Counts file: header line, then the genes in blocks (dense integer rows, as v1)
header <- paste(c("GENES", colnames(counts)), collapse = "\t")
writeLines(header, DATA_FILE)
n_genes <- nrow(counts)
block_starts <- seq(1, n_genes, by = GENES_PER_CHUNK)
for (block_start in block_starts) {
  block_end <- min(block_start + GENES_PER_CHUNK - 1, n_genes)
  block <- as.matrix(counts[block_start:block_end, , drop = FALSE])
  block_table <- data.table(GENES = rownames(block), block)
  fwrite(block_table, DATA_FILE, sep = "\t", append = TRUE, col.names = FALSE)
}
cat(sprintf("written %s (%d genes x %d cells) and %s\n", basename(DATA_FILE), n_genes, ncol(counts), basename(LABEL_FILE)))
