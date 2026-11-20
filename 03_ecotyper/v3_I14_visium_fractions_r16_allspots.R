# =============================================================================
# v3_I14_visium_fractions_r16_allspots.R — Visium spot fractions and 10x spot matrices for the EcoTyper Visium
# recovery, from ALL spots of the ROUND-16 CytoSPACE run, with the round-16 rule only
# (author 2026-11-10: "Visium (Fig 5b-e, S5b-d): use v3 rule, re run is possible").
#
# Why: v3_I4 started from the spot list of the v1 EcoTyper input (6,091 spots), which had already lost 81 spots by the
# v1 rule (v1 CytoSPACE cycling T cells > 0.5). 49 of those 81 do not meet the round-16 rule and 42 hold HRS cells.
# Here the spot list is every spot that the round-16 CytoSPACE run placed cells in (6,172), and the only filter is
# the round-16 one.
#   - source of the fractions: CytoSPACE fractional abundances per spot (16 round-16 types), Cytospace/output_r16
#   - filter: drop the spots whose round-16 Cycling fraction is > 0.5 (the 13-type model has no Cycling type)
#   - pools as the bulk 13 types: B_postGC = B_GC + B_Memory; NK = NK + ILC; Cycling dropped; rows renormalised to 1
#   - barcodes written with "-" changed to "." in the fractions (as v1 and v3_I4)
#   - spot counts: Cytospace/cHL_Visium_Input/ST_data_<id>.csv (all spots). The script stops unless these counts
#     equal the v1 10x matrix for every v1 spot, so the new matrix is the v1 matrix plus the missing spots.
#   - spot positions: the rows of the two v1 EcoTyper inputs (..._Input_Cyto and ..._Input_Ciber), which together
#     hold every spot; a spot present in both must have the same row.
# Run on Sherlock (R 4.2.0):  Rscript v3_I14_visium_fractions_r16_allspots.R   (submitted by
#   v3_I14_visium_recovery_r16_allspots.sh)
# Outputs: OAK Ecotyper/visium/2026-11-10_Visium_fractions_r16_allspots/<id>_fractions_r16_allspots.txt
#          OAK Ecotyper/visium/2026-11-10_Visium_fractions_r16_allspots/dropped_spots_cycling_allspots.tsv
#          OAK Ecotyper/visium/2026-11-10_Visium_Input_r16_allspots/<id>/{barcodes.tsv.gz, features.tsv.gz,
#              matrix.mtx.gz, tissue_positions_list.csv}
# Randomness: none.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

# Step 1. Folders, tumours and the spot filter
BRIAN_DIR       <- "${PROJECT_ROOT}"
CYTOSPACE_DIR   <- file.path(BRIAN_DIR, "Cytospace/output_r16")                      # round-16 CytoSPACE results (v3_I2)
ST_DATA_DIR     <- file.path(BRIAN_DIR, "Cytospace/cHL_Visium_Input")                # spot counts given to CytoSPACE
V1_CYTO_DIR     <- file.path(BRIAN_DIR, "Ecotyper/visium/2025-07-21_Visium_Input_Cyto")   # v1 input (6,091 spots)
V1_CIBER_DIR    <- file.path(BRIAN_DIR, "Ecotyper/visium/2025-07-21_Visium_Input_Ciber")  # other v1 input (positions only)
OUTPUT_DIR      <- file.path(BRIAN_DIR, "Ecotyper/visium/2026-11-10_Visium_fractions_r16_allspots")
NEW_INPUT_DIR   <- file.path(BRIAN_DIR, "Ecotyper/visium/2026-11-10_Visium_Input_r16_allspots")
TUMOURS         <- c("tumor1", "tumor2", "tumor3", "tumor4")
MAX_CYCLING     <- 0.5        # round-16 rule: drop spots with a Cycling fraction above 0.5 (author 2026-11-01)
if (dir.exists(OUTPUT_DIR) || dir.exists(NEW_INPUT_DIR)) {
  stop("an output folder exists; refusing to overwrite")
}
dir.create(OUTPUT_DIR, recursive = TRUE)
dir.create(NEW_INPUT_DIR, recursive = TRUE)

# A plain helper to read a gzipped text file line by line (fread() needs R.utils for .gz; Sherlock's R lacks it)
read_gz_lines <- function(path) {
  connection <- gzfile(path)
  lines <- readLines(connection)
  close(connection)
  lines
}

dropped_table <- list()

for (tumour in TUMOURS) {
  # Step 2. CytoSPACE fractions (16 types) of every spot that received cells
  cytospace_file <- file.path(CYTOSPACE_DIR, paste0(tumour, "_results"), "cytospace_results",
                              "fractional_abundances_by_spot.csv")
  cytospace <- fread(cytospace_file)
  n_cytospace_spots <- nrow(cytospace)

  # Step 3. Spot counts of all spots (genes in rows, spots in columns); spot order = the order of this file
  st_counts <- fread(file.path(ST_DATA_DIR, paste0("ST_data_", tumour, ".csv")))
  st_genes  <- st_counts[[1]]
  st_spots  <- colnames(st_counts)[-1]
  stopifnot(all(cytospace$SpotID %in% st_spots))
  spot_order <- st_spots[st_spots %in% cytospace$SpotID]       # CytoSPACE spots, in the ST_data order
  cytospace <- cytospace[match(spot_order, cytospace$SpotID)]

  # Step 4. Drop the spots with Cycling > 0.5 (round-16 rule) and record them
  is_cycling_spot <- cytospace$Cycling > MAX_CYCLING
  dropped_table[[tumour]] <- data.table(tumour = tumour, SpotID = cytospace$SpotID[is_cycling_spot],
                                        Cycling = cytospace$Cycling[is_cycling_spot])
  cytospace <- cytospace[!is_cycling_spot]
  kept_spots <- cytospace$SpotID

  # Step 5. Pools of the 13 EcoTyper types; Cycling dropped
  cytospace[, B_postGC  := B_GC + B_Memory]
  cytospace[, NK_pooled := NK + ILC]
  pooled_or_dropped <- c("SpotID", "B_GC", "B_Memory", "NK", "ILC", "Cycling", "B_postGC", "NK_pooled")
  keep_types <- c(setdiff(names(cytospace), pooled_or_dropped), "B_postGC", "NK_pooled")
  fractions <- as.matrix(cytospace[, ..keep_types])
  colnames(fractions)[colnames(fractions) == "NK_pooled"] <- "NK"
  stopifnot(ncol(fractions) == 13)

  # Step 6. Renormalise each spot to sum 1 (after the filter no spot can be all Cycling; checked anyway)
  spot_sums <- rowSums(fractions)
  stopifnot(all(spot_sums > 0))
  fractions <- fractions / spot_sums

  # Step 7. Write Mixture (barcode, "-" -> ".") + the 13 types
  output <- data.table(Mixture = gsub("-", ".", kept_spots), fractions)
  fwrite(output, file.path(OUTPUT_DIR, paste0(tumour, "_fractions_r16_allspots.txt")), sep = "\t")

  # Step 8. The v1 10x matrix, to check the ST_data counts against (v1 spots only)
  v1_folder <- file.path(V1_CYTO_DIR, tumour)
  v1_spots  <- read_gz_lines(file.path(v1_folder, "barcodes.tsv.gz"))
  v1_features <- read_gz_lines(file.path(v1_folder, "features.tsv.gz"))
  stopifnot(length(v1_features) == length(st_genes))                 # same number of genes, same order assumed below
  v1_matrix_lines <- read_gz_lines(file.path(v1_folder, "matrix.mtx.gz"))
  header_line <- v1_matrix_lines[1]
  v1_size <- as.integer(strsplit(v1_matrix_lines[2], " ")[[1]])      # genes, spots, entries
  stopifnot(v1_size[1] == length(st_genes), v1_size[2] == length(v1_spots))
  v1_entries <- fread(text = v1_matrix_lines[-(1:2)], header = FALSE, col.names = c("gene", "spot", "count"))
  v1_entries[, barcode := v1_spots[spot]]

  # Step 9. The ST_data counts as (gene number, barcode, count) entries, non-zero only
  st_matrix <- as.matrix(st_counts[, -1])
  rm(st_counts)
  nonzero <- which(st_matrix != 0, arr.ind = TRUE)
  st_entries <- data.table(gene = nonzero[, 1], barcode = st_spots[nonzero[, 2]], count = st_matrix[nonzero])
  rm(st_matrix, nonzero)

  # Step 10. Check: for the v1 spots, the ST_data entries are exactly the v1 matrix entries
  st_on_v1 <- st_entries[barcode %in% v1_spots]
  setorder(st_on_v1, barcode, gene)
  v1_sorted <- v1_entries[, .(gene, barcode, count)]
  setorder(v1_sorted, barcode, gene)
  stopifnot(nrow(st_on_v1) == nrow(v1_sorted))
  stopifnot(all(st_on_v1$gene == v1_sorted$gene), all(st_on_v1$barcode == v1_sorted$barcode),
            all(st_on_v1$count == v1_sorted$count))
  # the gene names of ST_data agree with the gene names of the v1 feature file (second column), where both are given
  v1_gene_names <- vapply(strsplit(v1_features, "\t"), function(fields) fields[2], character(1))
  n_same_gene_name <- sum(st_genes == v1_gene_names)
  cat(sprintf("%s: ST_data counts equal the v1 matrix on all %d v1 spots (%d entries); gene names equal in %d of %d rows\n",
              tumour, length(v1_spots), nrow(v1_sorted), n_same_gene_name, length(st_genes)))

  # Step 11. New 10x folder: features (the v1 file), barcodes, matrix and positions of the kept spots
  tumour_folder <- file.path(NEW_INPUT_DIR, tumour)
  dir.create(tumour_folder)
  file.copy(file.path(v1_folder, "features.tsv.gz"), file.path(tumour_folder, "features.tsv.gz"))

  #   barcodes: kept spots, in the ST_data order
  barcode_connection <- gzfile(file.path(tumour_folder, "barcodes.tsv.gz"), "w")
  writeLines(kept_spots, barcode_connection)
  close(barcode_connection)

  #   matrix: entries of the kept spots, spot column numbered 1..n in the barcode order (Matrix Market coordinates)
  new_entries <- st_entries[barcode %in% kept_spots]
  new_entries[, spot := match(barcode, kept_spots)]
  setorder(new_entries, spot, gene)
  plain_matrix_file <- file.path(tumour_folder, "matrix.mtx")
  writeLines(c(header_line, paste(length(st_genes), length(kept_spots), nrow(new_entries))), plain_matrix_file)
  fwrite(new_entries[, .(gene, spot, count)], plain_matrix_file, sep = " ", col.names = FALSE, append = TRUE)
  gzip_status <- system2("gzip", plain_matrix_file)
  stopifnot(gzip_status == 0)

  #   tissue positions: rows of the two v1 inputs together (no header, as v1); a spot in both must agree
  positions_cyto  <- fread(file.path(V1_CYTO_DIR, tumour, "tissue_positions_list.csv"), header = FALSE)
  positions_ciber <- fread(file.path(V1_CIBER_DIR, tumour, "tissue_positions_list.csv"), header = FALSE)
  in_both <- intersect(positions_cyto$V1, positions_ciber$V1)
  rows_cyto  <- as.data.frame(positions_cyto[match(in_both, V1)])
  rows_ciber <- as.data.frame(positions_ciber[match(in_both, V1)])
  stopifnot(isTRUE(all.equal(rows_cyto, rows_ciber, check.attributes = FALSE)))
  positions <- rbind(positions_cyto, positions_ciber[!V1 %in% positions_cyto$V1])
  stopifnot(all(kept_spots %in% positions$V1))
  positions <- positions[match(kept_spots, V1)]
  fwrite(positions, file.path(tumour_folder, "tissue_positions_list.csv"), col.names = FALSE)

  n_new_spots <- sum(!kept_spots %in% v1_spots)
  cat(sprintf("%s: %d CytoSPACE spots, %d dropped (Cycling > %.1f), %d kept (%d of them not in the v1 input); %d types (%s)\n",
              tumour, n_cytospace_spots, sum(is_cycling_spot), MAX_CYCLING, length(kept_spots), n_new_spots,
              ncol(fractions), paste(colnames(fractions), collapse = " ")))
  rm(st_entries, new_entries, v1_entries, v1_sorted, st_on_v1)
}

# Step 12. Table of the dropped spots, for the author
fwrite(rbindlist(dropped_table), file.path(OUTPUT_DIR, "dropped_spots_cycling_allspots.tsv"), sep = "\t")
cat("DONE\n")
