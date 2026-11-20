# =============================================================================
# v2_E1_steidl_symbols.R — plan A5 for the Steidl 2010 microarray (130 tumours, log2): map its gene symbols to GENCODE v47.
#
# Rule
#   - A symbol already in v47 is kept.
#   - Otherwise limma::alias2Symbol (org.Hs.eg.db, current HGNC) gives the current symbol. It is used only if it is
#     unique, in v47, and not already a row of the array (so no rows are summed or duplicated).
#
# Input:  v1 key p("steidl_microarray") = data/01_input/Steidl_microarray.txt
#         (identical to OAK Cibersort/bulk_input_2023to2025/051923_nejm_2010_microarray.txt)
# Output: data_v2/02_intermediate/steidl/2026-10-22_Steidl_microarray_v47.txt
#         data_v2/02_intermediate/steidl/2026-10-22_Steidl_microarray_v47_symbol_map.tsv
# Runs locally:  Rscript v2_E1_steidl_symbols.R
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(limma); library(org.Hs.eg.db) })
MANUSCRIPT_DIR <- "${BOX_ROOT}/paper/Manuscript"

# Step 1. Read the array (first column renamed to "gene") and the v47 gene symbols
array_data <- fread(file.path(MANUSCRIPT_DIR, "data/01_input/Steidl_microarray.txt"))
setnames(array_data, 1, "gene")
v47_symbols <- unique(fread(file.path(MANUSCRIPT_DIR, "data_v2/02_intermediate/ref/gene_table_v47.tsv"))$symbol)

# Step 2. Which array symbols are not in v47
old_symbols <- array_data$gene
new_symbols <- old_symbols
not_in_v47  <- !old_symbols %in% v47_symbols
symbols_to_map <- old_symbols[not_in_v47]

# Step 3. Current HGNC symbol of each of them (NA when alias2Symbol gives none or more than one)
candidates <- character(length(symbols_to_map))
for (k in seq_along(symbols_to_map)) {
  current <- suppressWarnings(alias2Symbol(symbols_to_map[k], species = "Hs"))
  if (length(current) == 1) {
    candidates[k] <- current
  } else {
    candidates[k] <- NA_character_
  }
}

# Step 4. Accept a candidate only if it is: found, in v47, not already an array symbol,
#   and proposed for exactly one array row (candidates that appear twice are both rejected)
accepted <- !is.na(candidates) &
            candidates %in% v47_symbols &
            !candidates %in% old_symbols &
            !duplicated(candidates) &
            !duplicated(candidates, fromLast = TRUE)
rows_to_rename <- which(not_in_v47)[accepted]
new_symbols[rows_to_rename] <- candidates[accepted]

# Step 5. Symbol map: in_v47 (kept), renamed, or not_in_v47 (kept, but not a v47 symbol)
status <- fifelse(!not_in_v47, "in_v47", fifelse(new_symbols != old_symbols, "renamed", "not_in_v47"))
symbol_map <- data.table(array_symbol = old_symbols, v47_symbol = new_symbols, status = status)

# Step 6. Write the renamed array and the map
array_data[, gene := new_symbols]
out_dir <- file.path(MANUSCRIPT_DIR, "data_v2/02_intermediate/steidl")
dir.create(out_dir, showWarnings = FALSE)
fwrite(array_data, file.path(MANUSCRIPT_DIR, "data_v2/02_intermediate/steidl/2026-10-22_Steidl_microarray_v47.txt"), sep = "\t")
fwrite(symbol_map, file.path(MANUSCRIPT_DIR, "data_v2/02_intermediate/steidl/2026-10-22_Steidl_microarray_v47_symbol_map.tsv"), sep = "\t")
print(symbol_map[, .N, by = status])

# Step 7. Check: how many signature-matrix genes and EcoTyper state genes are in the array, before -> after renaming
signature_genes <- fread(file.path(MANUSCRIPT_DIR, "data_v2/02_intermediate/stepD/merged/sigmtx.txt"))$NAME
gene_info_files <- list.files(file.path(MANUSCRIPT_DIR, "data_v2/02_intermediate/ecotyper_final"), "gene_info.txt",
                              recursive = TRUE, full.names = TRUE)
state_genes <- c()
for (gene_info_file in gene_info_files) {
  state_genes <- c(state_genes, fread(gene_info_file)$Gene)
}
state_genes <- unique(state_genes)
cat(sprintf("signature genes in array: %d -> %d of %d; EcoTyper state genes: %d -> %d of %d\n",
            sum(signature_genes %in% old_symbols), sum(signature_genes %in% new_symbols), length(signature_genes),
            sum(state_genes %in% old_symbols), sum(state_genes %in% new_symbols), length(state_genes)))
