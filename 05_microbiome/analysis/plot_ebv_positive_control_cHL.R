#!/usr/bin/env Rscript
# =============================================================================
# plot_ebv_positive_control_cHL.R
# Derived from the EBV positive-control panel in diff_abundance.R, but:
#   (1) NLPHL removed  -> compares EBV+ cHL vs EBV- cHL only
#   (2) adds a statistical test on EBV RPM (two-sided Wilcoxon rank-sum,
#       the study's primary zero-inclusive DA test) with a significance bracket.
# Output: analysis/ebv_concordance/EBV_positive_control_cHL.pdf
# =============================================================================
suppressMessages({library(ggplot2); library(scales)})
set.seed(1)
dir.create("analysis/ebv_concordance", recursive=TRUE, showWarnings=FALSE)

## ---- load (identical to diff_abundance.R) ----
read_cv <- function(f, cohort){ d <- read.csv(f, stringsAsFactors=FALSE); d$cohort <- cohort; d }
cv <- rbind(read_cv("data/cHL_cross_validated_species.csv","cHL"),
            read_cv("data/NLPHL_cross_validated_species.csv","NLPHL"))
design_all <- read.csv("data/sample_design_clean.csv", stringsAsFactors=FALSE)
design <- design_all[design_all$keep=="TRUE", ]
cv <- cv[cv$sample %in% design$sample, ]
cv <- cv[!cv$is_host & !cv$is_blocklisted, ]
cv$taxid <- as.character(cv$taxid)

EBV_TAXIDS <- c("3050299","10376")          # HHV-4 species node + S1

## ---- taxid x sample RPM matrix (absent = 0) ----
taxa <- unique(cv[,c("taxid","species_label","superkingdom","genus")])
taxa <- taxa[!duplicated(taxa$taxid),]; rownames(taxa) <- taxa$taxid
samples <- design$sample
M <- matrix(0, nrow=nrow(taxa), ncol=length(samples), dimnames=list(taxa$taxid, samples))
M[cbind(cv$taxid, cv$sample)] <- cv$RPM
depth <- tapply(cv$total_reads, cv$sample, function(x) x[1])[samples]

## ---- EBV RPM per sample, cHL only, two groups ----
ebv_row <- which(taxa$taxid %in% EBV_TAXIDS)[1]
grp <- ifelse(design$cohort=="cHL" & design$ebv=="pos", "EBV+ cHL",
        ifelse(design$cohort=="cHL" & design$ebv=="neg", "EBV- cHL", NA))
ebvdf <- data.frame(group=grp, rpm=M[ebv_row, design$sample])
ebvdf <- ebvdf[!is.na(ebvdf$group), ]                       # drop NLPHL + EBV-unknown cHL
ebvdf$group <- factor(ebvdf$group, levels=c("EBV+ cHL","EBV- cHL"))
pc_e <- 0.1                                                 # fixed display pseudocount (RPM); zeros drawn at 0.1 on log axis

## ---- statistical test: zero-inclusive two-sided Wilcoxon rank-sum on RPM ----
wt <- wilcox.test(rpm ~ group, data=ebvdf)                  # tested on raw RPM (rank-based; log is monotone)
p  <- wt$p.value
p_lab <- if (p < 2.2e-16) "p < 2.2e-16" else
         if (p < 1e-3)    sprintf("p = %.1e", p) else sprintf("p = %.3f", p)
n_pos <- sum(ebvdf$group=="EBV+ cHL"); n_neg <- sum(ebvdf$group=="EBV- cHL")
med <- tapply(ebvdf$rpm, ebvdf$group, median)
cat(sprintf("Wilcoxon EBV+ cHL (n=%d, median RPM=%.3g) vs EBV- cHL (n=%d, median RPM=%.3g): W=%.0f, %s\n",
            n_pos, med["EBV+ cHL"], n_neg, med["EBV- cHL"], unname(wt$statistic), p_lab))

## ---- significance bracket geometry (data units, y is log10) ----
ytop  <- max(ebvdf$rpm + pc_e)
y_bar <- ytop * 2.4       # horizontal bar
y_tip <- ytop * 1.5       # bracket tick bottoms
y_txt <- ytop * 3.6       # p-value / n label
brk <- data.frame(x=c(1,1,2,2), y=c(y_tip,y_bar,y_bar,y_tip))

pE <- ggplot(ebvdf, aes(group, rpm+pc_e, fill=group)) +
  geom_boxplot(outlier.shape=NA, width=.6) +
  geom_jitter(width=.12, size=.7, alpha=.45) +
  geom_path(data=brk, aes(x,y), inherit.aes=FALSE, linewidth=.4) +
  annotate("text", x=1.5, y=y_txt,
           label=sprintf("Wilcoxon rank-sum, %s", p_lab), size=4, fontface="bold") +
  scale_y_log10(breaks=c(0.1,1,10,100,1000,10000),
                minor_breaks=as.vector(outer(2:9, 10^(-1:4))),  # 2:9 -> minor ticks only, no decade overlap
                labels=comma, expand=expansion(mult=c(.05,.12))) +
  annotation_logticks(sides="l", outside=TRUE) + coord_cartesian(clip="off") +
  scale_x_discrete(labels=c("EBV+ cHL"=sprintf("EBV+ cHL\n(n=%d)", n_pos),
                            "EBV- cHL"=sprintf("EBV- cHL\n(n=%d)", n_neg))) +
  scale_fill_manual(values=c("EBV+ cHL"="#E7298A","EBV- cHL"="#FBB4AE"), guide="none") +
  labs(title="EBV / HHV-4 RPM by EBV status (cHL only)",
       subtitle=sprintf("Positive control HHV-4 (taxid 3050299); median RPM %.3g vs %.3g",
                        med["EBV+ cHL"], med["EBV- cHL"]),
       x=NULL, y="EBV RPM + pseudocount (log10)") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))

ggsave("analysis/ebv_concordance/EBV_positive_control_cHL.pdf", pE, width=6.2, height=5.5, device=cairo_pdf)
cat("WROTE analysis/ebv_concordance/EBV_positive_control_cHL.pdf\n")
