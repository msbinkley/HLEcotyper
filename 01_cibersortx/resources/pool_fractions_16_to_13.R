#!/usr/bin/env Rscript
# =============================================================================
# pool_fractions_16_to_13.R
# Collapse the 16 CIBERSORTx cell-type fractions into the 13 pooled cell types used
# for EcoTyper cell-state discovery (Methods; mapping in resources/cell_type_mapping_16_to_13.tsv):
#   B_cell  = B_preGC + B_Memory
#   NK_cell = NK + ILC
#   T_MKI67 (proliferating T cells) is dropped
# then each sample is renormalised so the 13 fractions sum to 1.
#
# Usage: Rscript pool_fractions_16_to_13.R CIBERSORTxGEP_NA_Fractions-Adjusted.txt cHL_Fractions_Norm.txt
# The output is the "Cell type fractions" file referenced by 02_ecotyper/discovery_config.yaml.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: pool_fractions_16_to_13.R <Fractions-Adjusted.txt> <out.txt>")
fr <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE)
fr <- fr[, setdiff(colnames(fr), c("P-value", "Correlation", "RMSE"))]

pooled <- data.frame(
  Mixture         = fr$Mixture,
  HRS             = fr$HRS,
  T_CD8           = fr$T_CD8,
  Plasma          = fr$Plasma,
  DC_Conventional = fr$DC_Conventional,
  T_gd            = fr$T_gd,
  T_Treg          = fr$T_Treg,
  B_GC            = fr$B_GC,
  T_CD4           = fr$T_CD4,
  DC_Plasmacytoid = fr$DC_Plasmacytoid,
  Macrophage      = fr$Macrophage,
  CAF             = fr$CAF,
  B_cell          = fr$B_preGC + fr$B_Memory,
  NK_cell         = fr$NK + fr$ILC,
  check.names = FALSE)

num <- pooled[, -1]
pooled[, -1] <- num / rowSums(num)          # renormalise to 1 after dropping T_MKI67
write.table(pooled, args[2], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %s: %d samples x %d pooled cell types\n", args[2], nrow(pooled), ncol(pooled) - 1))
