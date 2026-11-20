#!/usr/bin/env Rscript
# =============================================================================
# 03_fit_LOD95.R — limit of detection (95 % detection probability) of HRS cells in cfRNA.
#
#   1. background threshold = mean + 3 SD of the deconvoluted HRS fraction in the
#      no-spike-in background mixtures
#   2. each spike-in mixture is "detected" if its HRS fraction exceeds the threshold
#   3. logistic regression  Detected ~ log10(spike-in fraction)
#   4. LOD95 = spike-in fraction at which the model predicts P(detected) = 0.95
#
# Inputs: CIBERSORTxGEP_NA_Fractions-Adjusted.txt (from 02_cibersortx_hires_spikein.sh)
#         <prefix>_design.txt (from 01_make_spikein_mixtures.R)
# Usage:  Rscript 03_fit_LOD95.R CIBERSORTxGEP_NA_Fractions-Adjusted.txt LOD95_design.txt out_prefix
# Output: <out_prefix>_LOD95.txt (threshold, LOD95) and <out_prefix>_LOD95.pdf (Fig. 6d)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(scales) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: 03_fit_LOD95.R <Fractions-Adjusted.txt> <design.txt> <out_prefix>")

fr  <- read.delim(args[1], check.names = FALSE) %>% select(Mixture, HRS)
des <- read.delim(args[2])
d   <- inner_join(des, fr, by = "Mixture")

## 1. background threshold
bg     <- d %>% filter(Set == "background")
cutoff <- mean(bg$HRS) + 3 * sd(bg$HRS)

## 2-3. detection calls + logistic model on log10(spike-in)
sp <- d %>% filter(Set == "spike", SpikeIn > 0) %>% mutate(Detected = as.integer(HRS > cutoff))
model <- glm(Detected ~ log10(SpikeIn), data = sp, family = "binomial")

## 4. LOD95 on the original fraction scale
logit_p95  <- log(0.95 / 0.05)
LOD95      <- 10^((logit_p95 - coef(model)[1]) / coef(model)[2])
cat(sprintf("background HRS threshold = %.4f%%   LOD95 = %.3f%% spike-in\n", 100 * cutoff, 100 * LOD95))
write.table(data.frame(background_threshold = cutoff, LOD95 = LOD95,
                       n_background = nrow(bg), n_spikein = nrow(sp)),
            paste0(args[3], "_LOD95.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

## figure: detection rate vs spike-in fraction with the fitted curve
p <- ggplot(sp, aes(SpikeIn, Detected)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = LOD95, linetype = "dotted", linewidth = 1) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE, colour = "blue") +
  stat_summary(fun.data = "mean_cl_boot", geom = "pointrange") +
  scale_x_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1), labels = percent_format(accuracy = 0.01)) +
  scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
  labs(title = "Limit of detection of HRS cells in cfRNA (logistic regression)",
       subtitle = sprintf("LOD95 = %s | background threshold = %s",
                          percent(LOD95, accuracy = 0.01), percent(cutoff, accuracy = 0.01)),
       x = "Spike-in fraction of HRS cells (log scale)", y = "Detection rate") +
  theme_minimal()
ggsave(paste0(args[3], "_LOD95.pdf"), p, width = 6, height = 4.5)
