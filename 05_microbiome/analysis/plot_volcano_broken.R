#!/usr/bin/env Rscript
# Volcano, 2x2, with a PROPER broken y-axis via ggbreak::scale_y_break (canonical break symbol,
# single coherent axis per panel). EBV is off-scale (-log10 FDR ~300) in #2/#1/#4 -> those panels
# break the axis (0-13 detail, then ~285-300 for EBV); #3 (EBV not extreme) is a normal axis.
# ggbreak objects don't patchwork-compose, so each panel is converted with ggplotify::as.ggplot()
# and assembled; one shared legend grob is extracted and placed at the bottom.
suppressMessages({library(ggplot2); library(ggbreak); library(ggrepel); library(ggplotify); library(patchwork); library(grid)})

cmp <- data.frame(
  key   = c("2_EBVpos_vs_EBVneg_cHL","1_cHL_vs_NLPHL","3_EBVneg_cHL_vs_NLPHL","4_EBVpos_cHL_vs_NLPHL"),
  title = c("EBV+ vs EBV- cHL  [internal]","cHL (all) vs NLPHL  [batch-confounded]",
            "EBV- cHL vs NLPHL  [batch-confounded]","EBV+ cHL vs NLPHL  [batch-confounded]"),
  broken= c(TRUE,TRUE,FALSE,TRUE), stringsAsFactors=FALSE)
cols <- c("EBV (HHV-4)"="#E7298A","contaminant genus"="#D95F02","FDR<0.05"="#2C7FB8","ns"="grey78")
read1 <- function(k){ r<-read.csv(sprintf("table/diff_abundance/%s_results.csv",k),stringsAsFactors=FALSE)
  r$y <- -log10(pmax(r$wilcox_padj,1e-300)); r$sig <- r$wilcox_padj<0.05
  r$class <- factor(ifelse(r$is_EBV=="TRUE","EBV (HHV-4)", ifelse(r$contaminant=="TRUE","contaminant genus",
              ifelse(r$sig,"FDR<0.05","ns"))), levels=names(cols)); r }
xlim <- {a<-do.call(rbind,lapply(cmp$key,read1)); range(a$log2FC,na.rm=TRUE)+c(-.5,.5)}

panel <- function(i){
  r <- read1(cmp$key[i]); brk <- cmp$broken[i]
  top <- r[r$sig & r$class=="FDR<0.05" & r$y<=13,]; top <- head(top[order(-top$y),],8)
  p <- ggplot(r, aes(log2FC,y)) +
    geom_hline(yintercept=-log10(0.05), linetype=2, colour="grey60") +
    geom_vline(xintercept=0, colour="grey85", linewidth=.3) +
    geom_point(aes(colour=class, size=class=="EBV (HHV-4)"), alpha=.75, show.legend=FALSE) +
    geom_text_repel(data=top, aes(label=species), size=2.3, max.overlaps=Inf, min.segment.length=0,
                    segment.color="grey60", ylim=c(-Inf, 12)) +   # keep labels in the lower segment (avoid ggbreak gap duplication)
    scale_colour_manual(values=cols, drop=FALSE, limits=names(cols)) +
    scale_size_manual(values=c(`FALSE`=1.2,`TRUE`=3.2)) +
    coord_cartesian(xlim=xlim) +
    labs(title=cmp$title[i], x="log2 fold-change (mean RPM)", y="-log10 BH-FDR") +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(face="bold", size=9.5), legend.position="none",
          axis.text=element_text(size=8))
  if(brk){
    # guard: the hidden gap (13,285) must contain NO non-EBV points, and EBV must be above 285,
    # else the break would silently hide real data or clip EBV (reviewer-flagged robustness check).
    gap <- r$y>13 & r$y<285 & r$class!="EBV (HHV-4)"
    if(any(gap)) warning(sprintf("[%s] %d non-EBV point(s) fall in the hidden break gap (13,285)!", cmp$key[i], sum(gap)))
    ev <- r$y[r$class=="EBV (HHV-4)"]
    if(length(ev) && any(ev<=285)) warning(sprintf("[%s] EBV y=%.1f is at/below break-top 285", cmp$key[i], max(ev)))
    # EBV at ~300: break the empty 13..285 gap. ggbreak renders the canonical break symbol.
    p <- p + geom_text(data=r[r$class=="EBV (HHV-4)",], aes(label="EBV"),
                       colour=cols["EBV (HHV-4)"], fontface="bold", size=2.6, vjust=-0.9) +
         scale_y_break(c(13, 285), scales=0.32, space=0.10, ticklabels=c(290,300)) +
         scale_y_continuous(breaks=c(0,3,6,9,12))
  } else {
    p <- p + scale_y_continuous(breaks=scales::pretty_breaks(6))
  }
  ggplotify::as.ggplot(print(p))     # render ggbreak -> composable grob
}
panels <- lapply(seq_len(nrow(cmp)), panel)

## shared legend grob (extracted from a native ggplot with the full scale)
lp <- ggplot(data.frame(x=1, class=factor(names(cols),levels=names(cols))), aes(x,x,colour=class)) +
  geom_point(size=2.6) + scale_colour_manual(values=cols, name=NULL) +
  theme(legend.position="bottom", legend.text=element_text(size=9))
lg <- ggplotGrob(lp)                                  # robust to ggplot2 guide-box naming (guide-box / -bottom / -right)
gi <- which(grepl("^guide-box", lg$layout$name)); legend <- NULL
for(i in gi){ g <- lg$grobs[[i]]; if(!inherits(g,"zeroGrob")){ legend <- g; break } }
if(is.null(legend)) legend <- lg$grobs[[gi[1]]]
pleg <- ggplotify::as.ggplot(legend)

fig <- ((panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) / pleg) +
  plot_layout(heights=c(1,1,0.12)) +
  plot_annotation(title="Differential microbial abundance — broken y-axis where EBV is off-scale (#2,#1,#4)",
    caption="Broken panels: y-axis breaks the empty 13..285 region; EBV at -log10 FDR ~300 sits above the break. #3 = normal axis (EBV not extreme). +log2FC = higher in first-named group | dashed = FDR 0.05 | vs-NLPHL = batch-confounded",
    theme=theme(plot.title=element_text(face="bold", size=12), plot.caption=element_text(size=7, colour="grey45", hjust=0)))

ggsave("analysis/diff_abundance/volcano_all_comparisons.pdf", fig, width=14, height=11, device=cairo_pdf)
cat("WROTE analysis/diff_abundance/volcano_all_comparisons.pdf (ggbreak broken axis)\n")
