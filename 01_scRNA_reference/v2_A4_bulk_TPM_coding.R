# =============================================================================
# v2_A4_bulk_TPM_coding.R — v2 bulk TPM for CIBERSORTx and EcoTyper (Sherlock).
#
# Author decision 2026-10-21: the bulk keeps ALL v47 protein-coding genes, not only the 17,823-gene scRNA universe.
#   Genes   = v47 protein_coding + IG_C_gene + TR_C_gene (the IG/TR constant genes are in the universe and can be
#             signature genes); no PAR_Y copies.
#   Source  = salmon `-l ISR --allowDovetail` ENSG TPM; the 187 tumours (the columns of the v1 bulk TPM).
#   Scaling = renormalised to 1e6 per tumour (as v1's coding-only TPM).
# Outputs (in <W>/ref):
#   bulk_TPM_v47_coding_187_ensg.tsv    ENSG rows
#   bulk_TPM_v47_coding_187_symbol.tsv  symbol rows for CIBERSORTx / EcoTyper; ENSGs sharing a symbol (e.g. the
#                                       X/Y pseudoautosomal copies) are summed   <- USE THIS ONE for CIBERSORTx and EcoTyper
#   bulk_sample_check.tsv               identity check vs v1: log1p Spearman on shared symbols (moved here from v2_A)
# Run: Rscript v2_A4_bulk_TPM_coding.R   (optional env var TPM_ENSG = path of the salmon ENSG TPM table)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
W <- "${PROJECT_ROOT}/cHL_v2"            # project folder on Sherlock
TPM_ENSG <- Sys.getenv("TPM_ENSG", "${LAB_ROOT}/sequencing/pipeline/runs/rnaseq/Binkley_cHL/cHL_salmon_ISR_dovetail/cHL_salmon_ISR_dovetail_gene_ensg_TPM.tsv")

# stable(): drop the version from a gene ID, e.g. ENSG00000123.4 -> ENSG00000123 (a "_PAR_Y" ending is kept)
stable <- function(x) sub("\\.[0-9]+(_PAR_Y)?$", "\\1", x)

# Step 1. The coding genes (from the v47 gene table written by v2_A_gene_universe.R)
v47 <- fread(file.path(W, "ref/gene_table_v47.tsv"))
coding <- v47[type %in% c("protein_coding", "IG_C_gene", "TR_C_gene") & !par_y]

# Step 2. The 187 tumours = the sample columns of the v1 bulk TPM (header only)
samples <- setdiff(names(fread(file.path(W, "input/v1_bulk_TPM_187.txt"), nrows = 0)), "Gene")

# Step 3. Read the salmon TPM; add the stable gene ID
tpm <- fread(TPM_ENSG)
setnames(tpm, 1, "gene")
tpm[, ensg := stable(gene)]
stopifnot(length(samples) == 187, all(samples %in% names(tpm)), all(coding$ensg %in% tpm$ensg))

# Step 4. Keep the coding genes and the 187 tumours; rows with the same stable ID are summed
#         (one row per gene, in order of first appearance)
coding_tpm <- tpm[ensg %in% coding$ensg, c("ensg", samples), with = FALSE]
coding_tpm <- coding_tpm[, lapply(.SD, sum), by = ensg, .SDcols = samples]

# Step 5. Renormalise each tumour to 1e6; keep the share of the original TPM that the coding genes carried
tpm_matrix <- as.matrix(coding_tpm[, -1])
share <- colSums(tpm_matrix) / 1e6
tpm_matrix <- sweep(tpm_matrix, 2, colSums(tpm_matrix), "/") * 1e6
fwrite(data.table(Gene = coding_tpm$ensg, tpm_matrix), file.path(W, "ref/bulk_TPM_v47_coding_187_ensg.tsv"), sep = "\t")

# Step 6. ENSG rows -> symbol rows (ENSGs that share a symbol are summed; rows stay in order of first appearance)
symbols <- coding$symbol[match(coding_tpm$ensg, coding$ensg)]
symbol_tpm <- rowsum(tpm_matrix, symbols, reorder = FALSE)
fwrite(data.table(Gene = rownames(symbol_tpm), symbol_tpm), file.path(W, "ref/bulk_TPM_v47_coding_187_symbol.tsv"), sep = "\t")

# Step 7. Coverage summary: how many universe genes and merged-reference genes are present
univ <- fread(file.path(W, "ref/gene_universe.tsv"))
reference_file <- file.path(W, "output/stepC/reference/cHL_v2_scRNA_reference_4000_merged.txt")
ref_rows <- fread(cmd = sprintf("cut -f1 '%s'", reference_file))[[1]][-1]      # gene column, header line dropped
n_shared_symbols <- sum(duplicated(symbols))
cat(sprintf("coding genes: %d ENSG -> %d symbols (%d symbols shared by >1 ENSG, summed); 187 tumours; columns sum to 1e6\n",
            nrow(symbol_tpm) + n_shared_symbols, nrow(symbol_tpm), n_shared_symbols))
cat(sprintf("coding share of the original TPM: median %.1f%% (range %.1f-%.1f)\n", 100 * median(share), 100 * min(share), 100 * max(share)))
cat(sprintf("universe genes present: %d/%d; merged-reference genes present: %d/%d\n",
            sum(univ$symbol %in% rownames(symbol_tpm)), nrow(univ), sum(ref_rows %in% rownames(symbol_tpm)), length(ref_rows)))

# Step 8. Identity check vs v1: log1p Spearman between each tumour's v2 and v1 profile on the shared symbols
v1_bulk <- fread(file.path(W, "input/v1_bulk_TPM_187.txt"))
common <- intersect(rownames(symbol_tpm), v1_bulk$Gene)
v2_log <- log1p(symbol_tpm[common, samples])
v1_log <- log1p(as.matrix(v1_bulk[match(common, Gene), ..samples]))
v2_vs_v1 <- cor(v2_log, v1_log, method = "spearman")       # rows = v2 tumours, columns = v1 tumours
self <- diag(v2_vs_v1)
check <- data.table(sample = samples, self_rho = round(self, 4), best_match_in_v1 = samples[apply(v2_vs_v1, 1, which.max)])
check[, swap_suspect := best_match_in_v1 != sample]        # TRUE = another v1 tumour matches better than its own
fwrite(check, file.path(W, "ref/bulk_sample_check.tsv"), sep = "\t")
cat(sprintf("v1 genes present: %d/%d; self rho vs v1 on %d shared genes: median %.3f (min %.3f); best match not self: %d\n",
            sum(v1_bulk$Gene %in% rownames(symbol_tpm)), nrow(v1_bulk), length(common), median(self), min(self), sum(check$swap_suspect)))
