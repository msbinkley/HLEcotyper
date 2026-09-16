#!/usr/bin/env Rscript
# =============================================================================
# 01_make_spikein_mixtures.R — in-silico HRS spike-in for the cfRNA limit of detection.
#
# Two sets of simulated cfRNA profiles (Methods, "LOD95 analysis of HRS cells in cfRNA"):
#   BACKGROUND : healthy-volunteer cfRNA (volunteers CD11-CD20), no spike-in, N_REP noisy
#                replicates each -> defines the background HRS-fraction threshold
#   SPIKE-IN   : healthy-volunteer cfRNA (volunteers CD01-CD10) mixed with HRS cell-line
#                (L-428, LL-100 panel) TPM at 16 fractions from 100 % down to 0.001 %
#                (plus 0), N_REP noisy replicates each
# Mixture = frac * HRS_TPM + (1 - frac) * cfRNA_TPM, plus Gaussian noise
# (sd = 0.1 % of the mean), floored at 0; every column renormalised to 1e6.
# Both sets are concatenated into one mixture matrix for CIBERSORTx HiRes
# (02_cibersortx_hires_spikein.sh); the resulting HRS fractions are fitted in 03_fit_LOD95.R.
#
# Inputs (assumed to exist):
#   CFRNA_TPM  cfRNA TPM, genes x samples, healthy volunteers named CD01..CD20 (column Gene)
#   HRS_TPM    bulk RNA-seq TPM of HRS cell lines (LL-100 panel), column Gene + L_428
# Usage: Rscript 01_make_spikein_mixtures.R cfRNA_TPM.txt LL100_TPM.txt out_prefix
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tibble) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: 01_make_spikein_mixtures.R <cfRNA_TPM> <HRS_TPM> <out_prefix>")
CFRNA_TPM <- args[1]; HRS_TPM <- args[2]; OUT <- args[3]
set.seed(1)

HRS_cell  <- "L_428"
spike_vol <- sprintf("CD%02d", 1:10)          # volunteers used for spike-in mixtures
bg_vol    <- sprintf("CD%02d", 11:20)         # volunteers used for background
N_REP     <- 10
spike_in_fractions <- c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001,
                        0.0005, 0.0002, 0.0001, 0.00005, 0.00002, 0.00001, 0)

renorm <- function(df) { for (j in seq_along(df)) df[[j]] <- df[[j]] / sum(df[[j]]) * 1e6; df }
add_noise <- function(v) pmax(0, v + rnorm(length(v), 0, 0.001 * mean(v)))

cf  <- read.delim(CFRNA_TPM, check.names = FALSE) %>% select(Gene, all_of(c(spike_vol, bg_vol)))
hrs <- read.delim(HRS_TPM,   check.names = FALSE) %>% select(Gene, all_of(HRS_cell))
common <- intersect(cf$Gene, hrs$Gene)
cf  <- cf[match(common, cf$Gene), ];  hrs <- hrs[match(common, hrs$Gene), ]
cf[, -1] <- renorm(cf[, -1]);         hrs[, -1] <- renorm(hrs[, -1])

## ---- spike-in mixtures --------------------------------------------------------
spike <- list()
for (frac in spike_in_fractions) for (s in spike_vol) for (i in seq_len(N_REP)) {
  mix <- frac * hrs[[HRS_cell]] + (1 - frac) * cf[[s]]
  spike[[sprintf("Frac_%s_Sample%s_Rep_%d", frac, s, i)]] <- add_noise(mix)
}
spike <- renorm(as.data.frame(spike, check.names = FALSE))

## ---- background (no spike-in) --------------------------------------------------
bg <- list()
for (s in bg_vol) for (i in seq_len(N_REP)) bg[[sprintf("%s_Rep_%d", s, i)]] <- add_noise(cf[[s]])
bg <- renorm(as.data.frame(bg, check.names = FALSE))

## ---- combined CIBERSORTx mixture ---------------------------------------------
all <- renorm(cbind(spike, bg)) %>% add_column(Gene = common, .before = 1)
write.table(all, paste0(OUT, "_HRS_background_spikein_TPM.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- design table used by 03_fit_LOD95.R -------------------------------------
design <- rbind(
  data.frame(Mixture = colnames(spike), Sample = sub("^Frac_.*_Sample(CD\\d+)_Rep_\\d+$", "\\1", colnames(spike)),
             SpikeIn = as.numeric(sub("^Frac_([^_]+)_.*$", "\\1", colnames(spike))), Set = "spike"),
  data.frame(Mixture = colnames(bg), Sample = sub("_Rep_\\d+$", "", colnames(bg)), SpikeIn = 0, Set = "background"))
write.table(design, paste0(OUT, "_design.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote %d spike-in + %d background mixtures (%d genes)\n", ncol(spike), ncol(bg), length(common)))
