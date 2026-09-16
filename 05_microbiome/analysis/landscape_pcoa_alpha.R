#!/usr/bin/env Rscript
# =============================================================================
# landscape_pcoa_alpha.R -- community-level microbial LANDSCAPE, EBV+ vs EBV- cHL.
# Complements per-taxon DA + set enrichment with two community-level views:
#   (1) Beta-diversity ordination: Robust-Aitchison (rCLR) distance -> PCoA,
#       PERMANOVA (adonis2) for a location/centroid shift, guarded by betadisper
#       (PERMDISP; between-sample dispersion) so a spread difference isn't misread
#       as a location shift.
#   (2) Alpha diversity / dominance: Hill numbers q0/q1/q2, Wilcoxon EBV+ vs EBV-.
# Run TWICE: EBV-IN (whole community) and EBV-OUT (HHV-4 removed). If separation
# attenuates without HHV-4, the community-level separation is predominantly
# HHV-4-dependent (NOT proof the rest of the community is equivalent; p=0.076 is a
# conservative test since the larger EBV- arm is the more dispersed one, A&W 2013).
# cHL ONLY, grouped by clinical EBER status; same cohort/batch => NOT batch-confounded.
# Universe = DA prevalence set.
# =============================================================================
suppressMessages({library(vegan); library(ggplot2); library(scales); library(patchwork)})
set.seed(1)
dir.create("analysis/landscape", recursive=TRUE, showWarnings=FALSE)
dir.create("table/landscape", recursive=TRUE, showWarnings=FALSE)
PERM <- 9999
EBV_TAXIDS <- c("3050299","10376")
COLS <- c("EBV+ cHL"="#E7298A","EBV- cHL"="#FBB4AE")

## ---- load cHL only, build RPM matrix (absent=0) ----
cv <- read.csv("data/cHL_cross_validated_species.csv", stringsAsFactors=FALSE)
des <- read.csv("data/sample_design_clean.csv", stringsAsFactors=FALSE)
des <- des[des$keep=="TRUE" & des$cohort=="cHL" & des$ebv %in% c("pos","neg"), ]
cv  <- cv[cv$sample %in% des$sample & !cv$is_host & !cv$is_blocklisted, ]
cv$taxid <- as.character(cv$taxid)
taxa <- cv[!duplicated(cv$taxid), c("taxid","species_label","genus")]; rownames(taxa) <- taxa$taxid
samples <- des$sample
group <- factor(ifelse(des$ebv=="pos","EBV+ cHL","EBV- cHL"), levels=c("EBV+ cHL","EBV- cHL"))
names(group) <- samples
M <- matrix(0, nrow(taxa), length(samples), dimnames=list(taxa$taxid, samples))
M[cbind(cv$taxid, cv$sample)] <- cv$RPM
P <- M > 0

## ---- feature universe = DA prevalence set (>=10% in either arm) ----
A <- samples[group=="EBV+ cHL"]; B <- samples[group=="EBV- cHL"]
prevA <- rowMeans(P[,A,drop=FALSE]); prevB <- rowMeans(P[,B,drop=FALSE])
univ  <- rownames(M)[prevA>=0.10 | prevB>=0.10]
feat_in  <- union(univ, EBV_TAXIDS[EBV_TAXIDS %in% rownames(M)])   # EBV forced in
feat_out <- setdiff(univ, EBV_TAXIDS)                              # EBV removed
cat(sprintf("cHL samples: EBV+ %d, EBV- %d | features: universe=%d, EBV-in=%d, EBV-out=%d\n",
            length(A), length(B), length(univ), length(feat_in), length(feat_out)))

## ---- one community-level pass (beta + alpha) for a given feature set ----
run_landscape <- function(feats, tag){
  Msub <- M[feats, , drop=FALSE]
  keepS <- colSums(Msub) > 0                       # drop samples empty in this feature set
  dropped <- sum(!keepS)
  Msub <- Msub[, keepS, drop=FALSE]; g <- droplevels(group[keepS])
  Xs <- t(Msub)                                    # vegan wants samples x features
  # --- beta: robust-Aitchison (rCLR) distance -> PCoA + PERMANOVA + dispersion ---
  d  <- vegdist(Xs, method="robust.aitchison")
  cm <- cmdscale(d, k=2, eig=TRUE)
  ve <- round(100 * cm$eig[1:2] / sum(cm$eig[cm$eig>0]), 1)
  set.seed(1); ad <- adonis2(d ~ g, permutations=PERM, by="terms")
  R2 <- ad$R2[1]; pP <- ad$`Pr(>F)`[1]
  # PERMDISP: type="median" (robust default), bias.adjust=TRUE for unequal n (75 vs 92; Anderson et al. 2006)
  set.seed(1); bd <- betadisper(d, g, type="median", bias.adjust=TRUE)
  pt <- permutest(bd, permutations=PERM); pD <- pt$tab$`Pr(>F)`[1]
  disp <- tapply(bd$distances, g, mean)                 # mean distance to spatial median, per group
  ord <- data.frame(sample=rownames(Xs), PC1=cm$points[,1], PC2=cm$points[,2], group=g)
  # --- alpha: Hill q0 (richness), q1 (exp-Shannon), q2 (invSimpson) ---
  alp <- data.frame(sample=rownames(Xs), group=g,
                    q0=specnumber(Xs),
                    q1=exp(diversity(Xs, index="shannon")),
                    q2=diversity(Xs, index="invsimpson"))
  a_stats <- do.call(rbind, lapply(c("q0","q1","q2"), function(q){
    pv <- tryCatch(wilcox.test(alp[[q]] ~ alp$group)$p.value, error=function(e) NA)
    data.frame(dataset=tag, index=q, p=pv,
               med_pos=median(alp[[q]][alp$group=="EBV+ cHL"]),
               med_neg=median(alp[[q]][alp$group=="EBV- cHL"]))
  }))
  list(tag=tag, ord=ord, ve=ve, R2=R2, pP=pP, pD=pD, dropped=dropped,
       disp_pos=disp["EBV+ cHL"], disp_neg=disp["EBV- cHL"],
       nfeat=length(feats), alpha=alp, a_stats=a_stats,
       n_pos=sum(g=="EBV+ cHL"), n_neg=sum(g=="EBV- cHL"))
}
res_in  <- run_landscape(feat_in,  "EBV-in (whole community)")
res_out <- run_landscape(feat_out, "EBV-out (HHV-4 removed)")

## ---- PERMANOVA / dispersion summary table ----
beta_tab <- do.call(rbind, lapply(list(res_in,res_out), function(r) data.frame(
  dataset=r$tag, n_features=r$nfeat, n_pos=r$n_pos, n_neg=r$n_neg, samples_dropped=r$dropped,
  PC1_pct=r$ve[1], PC2_pct=r$ve[2], permanova_R2=round(r$R2,4), permanova_p=r$pP,
  betadisper_type="median", betadisper_bias_adjust=TRUE, betadisper_p=r$pD,
  disp_mean_EBVpos=round(r$disp_pos,3), disp_mean_EBVneg=round(r$disp_neg,3),
  more_dispersed=ifelse(r$disp_neg>r$disp_pos,"EBV-","EBV+"))))
write.csv(beta_tab, "table/landscape/landscape_beta_stats.csv", row.names=FALSE)
alpha_tab <- rbind(res_in$a_stats, res_out$a_stats); alpha_tab$p_BH <- p.adjust(alpha_tab$p,"BH")
write.csv(alpha_tab, "table/landscape/landscape_alpha_stats.csv", row.names=FALSE)

## ---- PCoA panels (EBV-in | EBV-out) ----
pcoa_panel <- function(r){
  ggplot(r$ord, aes(PC1, PC2, colour=group, fill=group)) +
    stat_ellipse(type="norm", level=0.95, geom="polygon", alpha=.12, colour=NA) +
    stat_ellipse(type="norm", level=0.95, linewidth=.5) +
    geom_point(size=1.7, alpha=.8) +
    scale_colour_manual(values=COLS, name=NULL) + scale_fill_manual(values=COLS, guide="none") +
    labs(title=r$tag,
         subtitle=sprintf("PERMANOVA R²=%.3f, p=%s  |  dispersion p=%s",
                          r$R2, format.pval(r$pP,eps=1e-4), format.pval(r$pD,eps=1e-4)),
         x=sprintf("PCo1 (%.1f%%)", r$ve[1]), y=sprintf("PCo2 (%.1f%%)", r$ve[2])) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"), legend.position="bottom")
}
pg <- (pcoa_panel(res_in) | pcoa_panel(res_out)) +
  plot_annotation(title="Community landscape: Robust-Aitchison PCoA, EBV+ vs EBV- cHL",
                  subtitle="Same cohort/batch (not confounded). Separation attenuates when HHV-4 is removed -> predominantly HHV-4-dependent (p=0.076 is not proof of equivalence).",
                  theme=theme(plot.title=element_text(face="bold", size=13))) +
  plot_layout(guides="collect") & theme(legend.position="bottom")
ggsave("analysis/landscape/landscape_pcoa.pdf", pg, width=11, height=6, device=cairo_pdf)

## ---- alpha-diversity companion (q0/q1/q2 x EBV-in/EBV-out) ----
# per-panel free scales via facet_wrap on an ordered combined factor (facet_grid ties
# y across columns, which would squish q1/q2 under q0's much larger richness scale).
idx_lab   <- c(q0="q0 richness", q1="q1 exp-Shannon", q2="q2 inverse-Simpson")
dshort    <- c("EBV-in (whole community)"="EBV-in", "EBV-out (HHV-4 removed)"="EBV-out")
panel_lvl <- as.vector(t(outer(names(dshort), names(idx_lab),
                function(d,q) paste0(dshort[d], " · ", idx_lab[q]))))  # EBV-in row, then EBV-out row
mkpanel <- function(dataset, index) factor(paste0(dshort[dataset]," · ", idx_lab[index]), levels=panel_lvl)
adf <- rbind(res_in$alpha,  res_out$alpha)
adf$dataset <- rep(c(res_in$tag, res_out$tag), c(nrow(res_in$alpha), nrow(res_out$alpha)))
adf_long <- do.call(rbind, lapply(c("q0","q1","q2"), function(q)
  data.frame(group=adf$group, value=adf[[q]], panel=mkpanel(adf$dataset, q))))
# p-value labels per panel
plab <- alpha_tab
plab$panel <- mkpanel(plab$dataset, plab$index)
plab$lab <- ifelse(plab$p_BH < 1e-4, "BH p<1e-04", sprintf("BH p=%.2g", plab$p_BH))
plab <- merge(plab, aggregate(value~panel, adf_long, max), by="panel")
pa <- ggplot(adf_long, aes(group, value, fill=group)) +
  geom_boxplot(outlier.shape=NA, width=.6) + geom_jitter(width=.12, size=.5, alpha=.35) +
  geom_text(data=plab, aes(x=1.5, y=value*1.10, label=lab), inherit.aes=FALSE, size=3.1) +
  facet_wrap(~panel, ncol=3, scales="free_y") +
  scale_y_continuous(expand=expansion(mult=c(.05,.15))) +
  scale_fill_manual(values=COLS, guide="none") +
  labs(title="Alpha diversity / dominance: EBV+ vs EBV- cHL",
       subtitle="Hill numbers, prevalence-filtered community. EBV-in q1/q2 collapse = HHV-4 dominance; EBV-out (background) shows no significant difference.",
       x=NULL, y="Hill diversity (effective # taxa)") +
  theme_bw(base_size=11) + theme(plot.title=element_text(face="bold"), strip.background=element_rect(fill="grey92"))
ggsave("analysis/landscape/landscape_alpha.pdf", pa, width=9, height=6.5, device=cairo_pdf)

## ---- console summary ----
cat("\n================ COMMUNITY LANDSCAPE (EBV+ vs EBV- cHL) ================\n")
print(beta_tab, row.names=FALSE)
cat("\n-- alpha (Wilcoxon, BH across 6 tests) --\n"); print(alpha_tab, row.names=FALSE, digits=3)
cat("\nWROTE analysis/landscape/{landscape_pcoa,landscape_alpha}.pdf and table/landscape/*.csv\n")
