# =============================================================================
# v2_D_pool13.R — CIBERSORTx 16-type fractions -> the 13 EcoTyper types (plan D1; author 2026-10-24).
#   B_GC + B_Memory -> B_postGC   (author 2026-10-24; B_GC and B_Memory signatures correlate r 0.86)
#   NK + ILC        -> NK         (author 2026-10-24; the pooled type keeps the name NK)
#   B_preGC kept as its own type
#   Cycling dropped (replaces v1's T_MKI67)
#   then every sample is renormalised so its 13 fractions sum to 1.
#   (Before 2026-10-24 the pools were B_cell = B_preGC + B_Memory and NK_cell = NK + ILC, with B_GC kept separate.)
#
# Usage: Rscript v2_D_pool13.R <CIBERSORTx_Adjusted.txt> <out_fractions13.txt>
#   (called by the bulk jobs of v2_D1_cibersortx.sh)
# Output columns: Mixture, the 10 unchanged types (input order), B_postGC, NK.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

# Step 1. Read the CIBERSORTx fractions (one row per sample)
args <- commandArgs(TRUE)
input_file  <- args[1]
output_file <- args[2]
fractions16 <- fread(input_file)

# Step 2. The cell-type columns = every column except Mixture and the three CIBERSORTx statistics.
#   Stop if there are not exactly 16, or if a type to pool / drop is missing.
cell_types <- setdiff(names(fractions16), c("Mixture", "P-value", "Correlation", "RMSE"))
stopifnot(length(cell_types) == 16,
          all(c("B_preGC", "B_GC", "B_Memory", "NK", "ILC", "Cycling") %in% cell_types))

# Step 3. Pooled columns (made as new columns first, so the 16-type NK column is not overwritten before use)
fractions16[, B_postGC  := B_GC + B_Memory]
fractions16[, NK_pooled := NK + ILC]

# Step 4. Keep the 10 untouched types (in input order) + the two pooled types; Cycling is dropped
pooled_or_dropped <- c("B_GC", "B_Memory", "NK", "ILC", "Cycling")
keep_types <- c(setdiff(cell_types, pooled_or_dropped), "B_postGC", "NK_pooled")
fractions13 <- as.matrix(fractions16[, ..keep_types])
colnames(fractions13)[colnames(fractions13) == "NK_pooled"] <- "NK"   # the pooled NK + ILC column is named NK
cycling_fraction <- fractions16$Cycling                 # kept only for the summary line below

# Step 5. Renormalise: divide each sample's fractions by their sum, so each row sums to 1
fractions13 <- fractions13 / rowSums(fractions13)

# Step 6. Write Mixture + the 13 types, tab-separated
fwrite(data.table(Mixture = fractions16$Mixture, fractions13), output_file, sep = "\t")
cat(sprintf("%d samples; 13 types; Cycling dropped: median %.4f, max %.4f of the 16-type fractions\n",
            nrow(fractions13), median(cycling_fraction), max(cycling_fraction)))
