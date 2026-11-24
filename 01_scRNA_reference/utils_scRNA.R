# =============================================================================
# utils_scRNA.R — helper functions shared by the scRNA preprocessing scripts (v1 helpers).
#
# Sourced by the v1 scripts 00_aoki_rebuild.R, 01_cHL_reference_pseudobulk.R and 02_external_pseudobulk.R
# (folder 01_scRNA/). The v2 scripts v2_C4_reference_pseudobulk.R, v2_C9_pseudobulk_merged.R and
# v2_C9b_pseudobulk_finallabels.R source the copy at <W>/script/utils_scRNA.R on Sherlock and use
# make_pseudobulk(), write_pseudobulk() and write_dense_tsv().
#
# Functions (names, arguments and return values must stay exactly as they are — other scripts call them):
#   work_dir(...)                                  folder for intermediate files
#   write_dense_tsv(m, file, header, ...)          write a sparse genes x cells matrix as a dense text table
#   rowsum_sparse(m, group)                        sparse version of base::rowsum()
#   make_pseudobulk(counts, labels, cell_type, ...) pseudo-bulk mixtures with known spike-in fractions
#   write_pseudobulk(df, file)                     write a pseudo-bulk table
#
# Why the matrices are written in row blocks:
#   All matrices are kept as sparse dgCMatrix and written to text 500 genes at a time, so the 2-4 GB dense
#   text files (CytoSPACE / EcoTyper inputs) never have to exist densely in RAM.
#
# Number formatting matches what the original notebook used for each file (argument `writer`):
#   writer = "fwrite" -> data.table::fwrite: same text as utils::write.table for INTEGER counts
#                        (incl. 1e+05 style), fast; used for raw UMI count matrices
#   writer = "base"   -> utils::write.table: needed for non-integer values (fwrite rounds the
#                        15th significant digit differently)
#   writer = "readr"  -> readr::write_tsv (shortest round-trip digits, no scientific notation)
# =============================================================================
suppressMessages({
  library(Matrix)
  library(data.table)
})


# -----------------------------------------------------------------------------
# work_dir(...): folder for intermediate files.
#   The folder is $HLE_WORK (default "../work") followed by the sub-folder names given in `...`.
#   It is created the first time it is asked for. Returns the folder path.
# -----------------------------------------------------------------------------
work_dir <- function(...) {
  folder <- file.path(Sys.getenv("HLE_WORK", "../work"), ...)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  folder
}


# -----------------------------------------------------------------------------
# write_dense_tsv(): write a genes x cells sparse matrix as a dense tab-delimited text file.
#   m            genes x cells matrix (sparse is fine)
#   file         output file
#   header       column names to write on the first line; the first entry is the name of the gene column,
#                so length(header) must be ncol(m) + 1
#   rownames_out gene names written in the first column (default: the matrix row names)
#   writer       "fwrite", "base" or "readr" (see the number-formatting note at the top of this file)
#   chunk        number of genes converted to dense and written at a time (default 500)
# Returns the file name (invisibly).
# -----------------------------------------------------------------------------
write_dense_tsv <- function(m, file, header, rownames_out = rownames(m),
                            writer = c("fwrite", "base", "readr"), chunk = 500) {
  writer <- match.arg(writer)
  stopifnot(length(header) == ncol(m) + 1, length(rownames_out) == nrow(m))

  # Step 1. Write the header line (this also creates / empties the file)
  connection <- file(file, open = "w")
  writeLines(paste(header, collapse = "\t"), connection)
  close(connection)

  # Step 2. Append the genes block by block: `chunk` genes per block
  block_starts <- seq(1, nrow(m), by = chunk)
  for (block_start in block_starts) {
    block_end  <- min(block_start + chunk - 1, nrow(m))
    block_rows <- block_start:block_end

    # make this block dense, with simple column names c1, c2, ... (column names are not written)
    dense_block <- as.matrix(m[block_rows, , drop = FALSE])
    block_table <- data.table::as.data.table(dense_block)
    data.table::setnames(block_table, paste0("c", seq_len(ncol(block_table))))

    # put the gene names in front as the first column
    block_table <- cbind(data.table::data.table(.g = rownames_out[block_rows]), block_table)

    # append the block with the chosen writer
    if (writer == "fwrite") {
      data.table::fwrite(block_table, file, sep = "\t", append = TRUE, col.names = FALSE, quote = FALSE)
    } else if (writer == "base") {
      utils::write.table(block_table, file, sep = "\t", append = TRUE, col.names = FALSE, row.names = FALSE,
                         quote = FALSE)
    } else {
      readr::write_tsv(block_table, file, append = TRUE, col_names = FALSE, na = "NA", escape = "none")
    }
    rm(dense_block, block_table)
  }
  invisible(gc())
  invisible(file)
}


# -----------------------------------------------------------------------------
# rowsum_sparse(): sparse equivalent of base::rowsum(as.matrix(m), group).
#   Adds up the rows of m that share the same group, and returns the summed rows in exactly the order
#   base::rowsum would (the sorted unique groups).
# -----------------------------------------------------------------------------
rowsum_sparse <- function(m, group) {
  # the group order base::rowsum() would use (run rowsum on a dummy matrix to get it)
  group_order <- rownames(rowsum(matrix(0, length(group), 1), group))
  group_factor <- factor(as.character(group), levels = group_order)

  # groups x rows indicator matrix: 1 where a row belongs to a group; multiplying sums the rows of each group
  indicator <- Matrix::fac2sparse(group_factor)
  summed <- indicator %*% m
  rownames(summed) <- group_order
  methods::as(summed, "CsparseMatrix")
}


# -----------------------------------------------------------------------------
# make_pseudobulk(): pseudo-bulk mixtures with known spike-in fractions of one cell type.
#   Ported from the notebook chunks "2025-06-12 PseudoBulk" and "Make pseudo bulk for Stewart and 10x and Aoki".
#
#   counts     genes x cells count matrix
#   labels     cell type of every cell (same order as the columns of counts)
#   cell_type  the cell type X to spike in
#   fractions  spike-in fractions (default 0, 0.02, ..., 0.40)
#   n_boot     number of bootstrap draws per fraction (default 10)
#
# How it works, for a cell type X with n_X cells among N cells:
#   - the other N - n_X cells ("the remainder") are summed once;
#   - for every fraction f, round(f * (N - n_X) / (1 - f)) X cells are drawn WITH replacement (bootstrap),
#     n_boot = 10 times with set.seed(1), set.seed(2), ..., set.seed(10); their summed counts are averaged
#     over the 10 draws and added to the remainder;
#   - column "sample_control" = all N cells summed;
#   - every column is then scaled to counts per million.
# Column names: X + the fraction printed as "%.3f" without "0." (X000, X020, ... X400).
# Returns a data frame: column "gene" followed by the pseudo-bulk columns.
#
# Speed-up vs. the notebook (results identical): rowSums(counts[, idx]) is computed as counts %*% tabulate(idx),
# which gives the same integer sums; the random-number calls are unchanged.
# -----------------------------------------------------------------------------
make_pseudobulk <- function(counts, labels, cell_type, fractions = seq(0, 0.4, by = 0.02),
                            n_boot = 10) {
  stopifnot(ncol(counts) == length(labels))

  # Step 1. Split the cells into the spike-in type and the remainder
  cluster_cell_coordinate <- which(labels == cell_type)               # column positions of the X cells
  remain <- setdiff(seq_len(ncol(counts)), cluster_cell_coordinate)   # column positions of all other cells
  num_remain_cells <- length(remain)
  counts_remain <- Matrix::rowSums(counts[, remain, drop = FALSE])    # remainder summed once

  # Step 2. Control column: all cells summed
  bulk <- data.frame(row.names = rownames(counts))
  bulk$sample_control <- Matrix::rowSums(counts)

  # Step 3. One column per spike-in fraction
  for (fraction in fractions) {
    # number of X cells needed so that X makes up `fraction` of the mixture
    num_of_cluster_cells <- round(fraction * num_remain_cells / (1 - fraction))

    # draw that many X cells with replacement, n_boot times (seed = draw number 1..n_boot)
    boot <- matrix(0, nrow(counts), n_boot)
    for (seed in seq_len(n_boot)) {
      set.seed(seed)
      bootstrap_cells <- sample(cluster_cell_coordinate, num_of_cluster_cells, replace = TRUE)
      times_drawn <- tabulate(bootstrap_cells, nbins = ncol(counts))   # how often each cell was drawn
      boot[, seed] <- as.vector(counts %*% times_drawn)                # = summed counts of the drawn cells
    }
    counts_cluster <- rowSums(boot) / n_boot                           # average over the draws

    # column name, e.g. "HRS020" for fraction 0.02 (spaces in the type name become "_")
    colname <- paste0(gsub(" ", "_", cell_type), gsub("0\\.", "", sprintf("%.3f", fraction)))
    bulk[, colname] <- rowSums(cbind(counts_remain, counts_cluster))
  }

  # Step 4. Counts -> counts per million (same arithmetic as the notebook)
  for (col in colnames(bulk)) {
    newsum <- sum(bulk[, col])
    bulk[, col] <- (1000000 / newsum) * bulk[, col]
  }
  cbind(gene = rownames(bulk), bulk)
}


# -----------------------------------------------------------------------------
# write_pseudobulk(): write a pseudo-bulk table exactly as the notebook did
#   (tab-separated, no quotes, with a header line, no row names).
# -----------------------------------------------------------------------------
write_pseudobulk <- function(df, file) {
  write.table(df, file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}
