#!/usr/bin/env Rscript
# =============================================================================
# plot_volcano_EBV_cHL.R -- single-comparison Wilcoxon volcano: EBV+ vs EBV- cHL,
# with EBV / HHV-4 EXCLUDED (we are not comparing EBV abundance here).
# Source of stats: table/diff_abundance/2_EBVpos_vs_EBVneg_cHL_results.csv (zero-inclusive
# Wilcoxon-on-RPM + BH, from diff_abundance.R). +log2FC = higher in EBV+ cHL.
# EBV is the off-scale point that forced the broken axis in the 4-panel figure;
# dropping it leaves max -log10(FDR)~3.8 -> a normal linear axis.
# =============================================================================
suppressMessages({library(ggplot2); library(ggrepel); library(scales)})
dir.create("analysis/diff_abundance", recursive=TRUE, showWarnings=FALSE)

r <- read.csv("table/diff_abundance/2_EBVpos_vs_EBVneg_cHL_results.csv", stringsAsFactors=FALSE)
istrue <- function(x) x %in% c(TRUE, "TRUE")
r <- r[!istrue(r$is_EBV), ]                     # drop EBV/HHV-4: not comparing EBV abundance
r$y   <- -log10(pmax(r$wilcox_padj, 1e-300))
r$sig <- r$wilcox_padj < 0.05
cols <- c("FDR<0.05"="#2C7FB8","ns"="grey78")            # color by significance only (no contaminant category)
r$class <- factor(ifelse(r$sig, "FDR<0.05", "ns"), levels=names(cols))

## label the top significant taxa (by FDR)
top <- r[r$sig, ]; top <- head(top[order(top$wilcox_padj), ], 12)
n_pos <- sum(r$sig & r$log2FC>0)     # higher in EBV+ cHL
n_neg <- sum(r$sig & r$log2FC<0)     # higher in EBV- cHL
xlim <- max(abs(r$log2FC), na.rm=TRUE) * c(-1,1) * 1.05

p <- ggplot(r, aes(log2FC, y)) +
  geom_hline(yintercept=-log10(0.05), linetype=2, colour="grey60") +
  geom_vline(xintercept=0, colour="grey85", linewidth=.3) +
  geom_point(aes(colour=class), size=1.7, alpha=.75) +
  geom_text_repel(data=top, aes(label=species), size=2.6, max.overlaps=Inf,
                  min.segment.length=0, segment.color="grey60") +
  annotate("text", x=xlim[2]*0.92, y=max(r$y)*0.02, label="higher in EBV+ cHL →",
           hjust=1, size=3, colour="grey40") +
  annotate("text", x=xlim[1]*0.92, y=max(r$y)*0.02, label="← higher in EBV- cHL",
           hjust=0, size=3, colour="grey40") +
  scale_colour_manual(values=cols, drop=FALSE, name=NULL) +
  scale_y_continuous(breaks=pretty_breaks(6), expand=expansion(mult=c(.02,.08))) +
  coord_cartesian(xlim=xlim) +
  labs(title="Differential microbial abundance: EBV+ vs EBV- cHL",
       subtitle=sprintf("Wilcoxon rank-sum on RPM (BH-FDR); EBV / HHV-4 excluded.\n%d taxa at FDR<0.05: %d higher in EBV+, %d higher in EBV-.",
                        n_pos+n_neg, n_pos, n_neg),
       x="log2 fold-change (mean RPM)", y="-log10 BH-FDR") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(face="bold"), legend.position="bottom")

ggsave("analysis/diff_abundance/volcano_EBVpos_vs_EBVneg_cHL.pdf", p, width=8.5, height=6.5, device=cairo_pdf)
cat(sprintf("WROTE analysis/diff_abundance/volcano_EBVpos_vs_EBVneg_cHL.pdf   [n=%d taxa, EBV excluded; %d sig at FDR<0.05]\n",
            nrow(r), n_pos+n_neg))
