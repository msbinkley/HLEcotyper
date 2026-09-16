#!/usr/bin/env Rscript
# =============================================================================
# diff_abundance.R -- differential tumor-microbiome abundance, Hodgkin lymphoma
# Method (codex-agreed, 2026-06-26): prevalence-filtered (>=10% in either arm),
# zero-inclusive Wilcoxon rank-sum on RPM + Fisher presence/absence companion;
# RPM (no CLR/rarefaction); pseudocount only for log2FC/plots; per-comparison BH.
# Comparisons 1/3/4 are CROSS-COHORT = batch-confounded -> DISCOVERY ONLY.
# EBV/HHV-4 = positive control; plant-fungi/kitome genera = negative controls.
# =============================================================================
suppressMessages({library(ggplot2); library(scales); library(patchwork)})
set.seed(1)
dir.create("table/diff_abundance", recursive=TRUE, showWarnings=FALSE); dir.create("analysis/ebv_concordance", recursive=TRUE, showWarnings=FALSE)

## ---- load ----
read_cv <- function(f, cohort){
  d <- read.csv(f, stringsAsFactors=FALSE)
  d$cohort <- cohort; d
}
cv <- rbind(read_cv("data/cHL_cross_validated_species.csv","cHL"),
            read_cv("data/NLPHL_cross_validated_species.csv","NLPHL"))
# Cleaned design: NLPHL-misclassified-in-cHL removed + cHL multi-sample patients deduped
# (keep highest-depth). NLPHL dedup still PENDING user input -> comparisons 1/3/4 provisional.
design_all <- read.csv("data/sample_design_clean.csv", stringsAsFactors=FALSE)
design <- design_all[design_all$keep=="TRUE", ]
cat(sprintf("Design: %d kept; %d excluded (dup/misclassified, see exclude_reason): cHL %d, NLPHL %d\n",
            nrow(design), sum(design_all$keep!="TRUE"),
            sum(design_all$keep!="TRUE" & design_all$cohort=="cHL"),
            sum(design_all$keep!="TRUE" & design_all$cohort=="NLPHL")))
cv <- cv[cv$sample %in% design$sample, ]      # drop excluded samples from the call table

## ---- feature universe: drop host + kitome-blocklisted ----
cv <- cv[!cv$is_host & !cv$is_blocklisted, ]
# one row per (sample, taxid); RPM is depth-normalized (reads / total_reads*1e6)
cv$taxid <- as.character(cv$taxid)

## ---- EBV + contaminant tags (taxid/label robust) ----
EBV_TAXIDS <- c("3050299","10376")          # species node in step-05 + S1 (HHV-4)
is_ebv_row <- cv$taxid %in% EBV_TAXIDS | grepl("humangamma4|gammaherpesvirus 4", cv$species_label, ignore.case=TRUE)
# known reagent/environmental contaminant genera (plant fungi, skin, kitome bacteria)
CONTAM_GENERA <- c("Malassezia","Pyrenophora","Zymoseptoria","Cercospora","Ascochyta",
  "Alternaria","Penicillium","Aspergillus","Cladosporium","Fusarium","Botrytis","Puccinia",
  "Ustilago","Sclerotinia","Magnaporthe","Colletotrichum","Verticillium","Psilocybe",
  "Ralstonia","Bradyrhizobium","Cutibacterium","Sphingomonas","Pseudomonas","Acidovorax",
  "Methylobacterium","Stenotrophomonas","Herbaspirillum","Delftia","Variovorax","Streptomyces")

## ---- build taxid x sample RPM matrix (absent = 0) ----
taxa <- unique(cv[,c("taxid","species_label","superkingdom","genus")])
taxa <- taxa[!duplicated(taxa$taxid),]
rownames(taxa) <- taxa$taxid
samples <- design$sample
M <- matrix(0, nrow=nrow(taxa), ncol=length(samples), dimnames=list(taxa$taxid, samples))
M[cbind(cv$taxid, cv$sample)] <- cv$RPM
# presence matrix
P <- M > 0
# per-sample sequencing depth (for pseudocount = half a one-read RPM)
depth <- tapply(cv$total_reads, cv$sample, function(x) x[1])
depth <- depth[samples]

## ---- comparison definitions ----
arms <- function(cond) design$sample[which(cond)]   # which() drops NA (EBV-unknown samples)
cmp <- list(
  list(id="2_EBVpos_vs_EBVneg_cHL", lab="EBV+ cHL vs EBV- cHL", confound=FALSE,
       A=arms(design$cohort=="cHL" & design$ebv=="pos"), Aname="EBV+ cHL",
       B=arms(design$cohort=="cHL" & design$ebv=="neg"), Bname="EBV- cHL"),
  list(id="3_EBVneg_cHL_vs_NLPHL", lab="EBV- cHL vs NLPHL", confound=TRUE,
       A=arms(design$cohort=="cHL" & design$ebv=="neg"), Aname="EBV- cHL",
       B=arms(design$cohort=="NLPHL"), Bname="NLPHL"),
  list(id="4_EBVpos_cHL_vs_NLPHL", lab="EBV+ cHL vs NLPHL", confound=TRUE,
       A=arms(design$cohort=="cHL" & design$ebv=="pos"), Aname="EBV+ cHL",
       B=arms(design$cohort=="NLPHL"), Bname="NLPHL"),
  list(id="1_cHL_vs_NLPHL", lab="cHL (all) vs NLPHL", confound=TRUE,
       A=arms(design$cohort=="cHL"), Aname="cHL",
       B=arms(design$cohort=="NLPHL"), Bname="NLPHL")
)

PREV_MIN <- 0.10
run_cmp <- function(cc){
  A <- cc$A; B <- cc$B
  pc <- 0.5 * 1e6 / median(depth[c(A,B)], na.rm=TRUE)   # half a one-read RPM
  prevA <- rowMeans(P[,A,drop=FALSE]); prevB <- rowMeans(P[,B,drop=FALSE])
  keep <- (prevA>=PREV_MIN | prevB>=PREV_MIN) | (taxa$taxid %in% EBV_TAXIDS)
  ti <- which(keep)
  res <- lapply(ti, function(i){
    a <- M[i,A]; b <- M[i,B]
    if (sum(a)==0 && sum(b)==0) return(NULL)
    wp <- tryCatch(wilcox.test(a,b)$p.value, error=function(e) NA_real_)
    # Fisher on presence/absence
    ft <- tryCatch(fisher.test(matrix(c(sum(a>0),sum(a==0),sum(b>0),sum(b==0)),2))$p.value,
                   error=function(e) NA_real_)
    data.frame(taxid=taxa$taxid[i], species=taxa$species_label[i],
               superkingdom=taxa$superkingdom[i], genus=taxa$genus[i],
               nA=length(A), nB=length(B), prevA=mean(a>0), prevB=mean(b>0),
               meanRPM_A=mean(a), meanRPM_B=mean(b), medRPM_A=median(a), medRPM_B=median(b),
               log2FC=log2((mean(a)+pc)/(mean(b)+pc)),
               wilcox_p=wp, fisher_p=ft, stringsAsFactors=FALSE)
  })
  res <- do.call(rbind, res)
  res$wilcox_padj <- p.adjust(res$wilcox_p, "BH")
  res$fisher_padj <- p.adjust(res$fisher_p, "BH")
  res$contaminant <- res$genus %in% CONTAM_GENERA
  res$is_EBV <- res$taxid %in% EBV_TAXIDS
  res <- res[order(res$wilcox_padj, -abs(res$log2FC)),]
  attr(res,"meta") <- cc; res
}

all_res <- lapply(cmp, run_cmp); names(all_res) <- sapply(cmp,function(x)x$id)

## ---- write tables ----
for(id in names(all_res)){
  write.csv(all_res[[id]], sprintf("table/diff_abundance/%s_results.csv", id), row.names=FALSE)
}

## ---- volcano per comparison ----
volcano <- function(res){
  cc <- attr(res,"meta")
  res$sig <- res$wilcox_padj < 0.05
  res$y <- -log10(pmax(res$wilcox_padj, 1e-300))
  res$col <- ifelse(res$is_EBV,"EBV (HHV-4)", ifelse(res$contaminant,"contaminant genus",
              ifelse(res$sig,"FDR<0.05","ns")))
  top <- head(res[res$sig & !res$contaminant,], 12)
  ttl <- sprintf("%s%s", cc$lab, if(cc$confound) "  [batch-confounded: discovery only]" else "  [internal]")
  ggplot(res, aes(log2FC, y)) +
    geom_hline(yintercept=-log10(0.05), linetype=2, colour="grey60") +
    geom_vline(xintercept=0, colour="grey80") +
    geom_point(aes(colour=col, size=is_EBV), alpha=.7) +
    ggrepel_or_text(top) +
    scale_colour_manual(values=c("EBV (HHV-4)"="#E7298A","contaminant genus"="#D95F02",
        "FDR<0.05"="#2C7FB8","ns"="grey75"), name=NULL) +
    scale_size_manual(values=c(`FALSE`=1.3,`TRUE`=3.2), guide="none") +
    labs(title=ttl, subtitle=sprintf("%s (n=%d) vs %s (n=%d) | +log2FC = higher in %s",
         cc$Aname,length(cc$A),cc$Bname,length(cc$B),cc$Aname),
         x="log2 fold-change (mean RPM)", y="-log10 BH-FDR") +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold", size=12))
}
# label helper: use ggrepel if present else geom_text
ggrepel_or_text <- function(top){
  if(nrow(top)==0) return(NULL)
  if(requireNamespace("ggrepel", quietly=TRUE))
    ggrepel::geom_text_repel(data=top, aes(label=species), size=2.6, max.overlaps=20, min.segment.length=0)
  else geom_text(data=top, aes(label=species), size=2.4, vjust=-.6, check_overlap=TRUE)
}

# NOTE: volcano_all_comparisons.pdf is generated by the dedicated plot_volcano_broken.R
# (per-panel broken y-axis, since EBV is off-scale). Do NOT write it here or it clobbers that
# version on re-run. The per-comparison volcano() above is retained only as a reference impl.

## ---- EBV positive-control panel (RPM across clinical groups) ----
ebv_row <- which(taxa$taxid %in% EBV_TAXIDS)[1]
grp <- factor(ifelse(design$cohort=="NLPHL","NLPHL",
              ifelse(design$ebv=="pos","EBV+ cHL", ifelse(design$ebv=="neg","EBV- cHL","cHL EBV-NA"))),
              levels=c("EBV+ cHL","EBV- cHL","cHL EBV-NA","NLPHL"))
ebvdf <- data.frame(group=grp, rpm=M[ebv_row, design$sample])
ebvdf <- ebvdf[!is.na(ebvdf$group) & ebvdf$group!="cHL EBV-NA",]; ebvdf$group <- droplevels(ebvdf$group)
pc_e <- 0.5*1e6/median(depth, na.rm=TRUE)
pE <- ggplot(ebvdf, aes(group, rpm+pc_e, fill=group)) +
  geom_boxplot(outlier.shape=NA, width=.6) +   # outliers shown via the jitter layer instead (avoids double-plotting)
  geom_jitter(width=.12, size=.6, alpha=.4) +
  scale_y_log10(labels=comma) +
  scale_fill_manual(values=c("EBV+ cHL"="#E7298A","EBV- cHL"="#FBB4AE","NLPHL"="#2C7FB8"), guide="none") +
  labs(title="EBV / HHV-4 positive control (taxid 3050299)",
       subtitle="EBV+ cHL >> EBV- cHL > NLPHL — supports EBV signal recovery + cHL metadata join (NLPHL EBV status is a caveat, not validated)",
       x=NULL, y="EBV RPM + pseudocount (log10)") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))
ggsave("analysis/ebv_concordance/EBV_positive_control.pdf", pE, width=7, height=5.5, device=cairo_pdf)

## ---- console summary ----
cat("\n================ DIFFERENTIAL ABUNDANCE SUMMARY ================\n")
for(id in names(all_res)){
  r <- all_res[[id]]; cc <- attr(r,"meta")
  sig <- sum(r$wilcox_padj<0.05, na.rm=TRUE); signc <- sum(r$wilcox_padj<0.05 & !r$contaminant, na.rm=TRUE)
  cat(sprintf("\n[%s] %s%s\n  tested=%d  FDR<0.05=%d  (non-contaminant=%d)\n",
      id, cc$lab, if(cc$confound)"  *DISCOVERY-ONLY*" else "", nrow(r), sig, signc))
  e <- r[r$is_EBV,]
  if(nrow(e)) cat(sprintf("  EBV: prevA=%.2f prevB=%.2f log2FC=%.2f padj=%.2g\n",
      e$prevA[1],e$prevB[1],e$log2FC[1],e$wilcox_padj[1]))
  top <- head(r[r$wilcox_padj<0.05 & !r$contaminant,c("species","log2FC","prevA","prevB","wilcox_padj")],6)
  if(nrow(top)) { cat("  top non-contaminant hits:\n"); print(top, row.names=FALSE) }
}
cat("\nWROTE table/diff_abundance/*_results.csv and analysis/ebv_concordance/*.pdf\n")
