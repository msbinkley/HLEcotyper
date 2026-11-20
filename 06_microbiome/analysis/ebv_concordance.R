#!/usr/bin/env Rscript
# =============================================================================
# ebv_concordance.R -- concordance: clinical EBER/IHC vs microbial-pipeline EBV.
# Method (codex-agreed): primary pipeline-EBV+ rule = EBV RPM >= 1 (RPM, not raw
# reads, to remove depth confound). Sensitivity rules: Youden-optimal RPM, RPM>=5,
# pipeline high_confidence flag, any-read. Clinical EBER = imperfect reference
# (discordants, NOT "false positives"). Metrics with Wilson CIs; kappa w/ bootstrap;
# AUC w/ DeLong CI; McNemar for paired rule comparison. Cleaned cHL only (n=199).
# =============================================================================
suppressMessages({library(ggplot2); library(scales); library(patchwork); library(pROC)})
set.seed(1); dir.create("analysis/ebv_concordance", recursive=TRUE, showWarnings=FALSE); dir.create("table/ebv_concordance", recursive=TRUE, showWarnings=FALSE)
d <- read.csv("data/ebv_concordance_input.csv", stringsAsFactors=FALSE)
gt <- d[d$clinical_ebv %in% c("pos","neg"), ]           # ground-truth subset
truth <- gt$clinical_ebv=="pos"; score <- gt$ebv_rpm
cat(sprintf("Ground-truth cHL: n=%d (pos=%d neg=%d); clinical-NA classified separately=%d\n",
            nrow(gt), sum(truth), sum(!truth), sum(is.na(d$clinical_ebv))))

## ---- ROC / AUC (DeLong) + Youden ----
roc1 <- roc(truth, score, quiet=TRUE, direction="<")
auc_ci <- ci.auc(roc1, method="delong")
yj <- coords(roc1, "best", best.method="youden", ret=c("threshold"))
thr_youden <- as.numeric(yj$threshold[1]); if(!is.finite(thr_youden)) thr_youden <- 1.7
cat(sprintf("AUC = %.3f (DeLong 95%% CI %.3f-%.3f); Youden threshold = %.2f RPM\n",
            as.numeric(auc(roc1)), auc_ci[1], auc_ci[3], thr_youden))

## ---- helpers ----
wilson <- function(x,n){ if(n==0) return(c(NA,NA)); z<-1.96; p<-x/n; d<-1+z^2/n
  c((p+z^2/(2*n)-z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d,(p+z^2/(2*n)+z*sqrt(p*(1-p)/n+z^2/(4*n^2)))/d) }
ckappa <- function(tr,ca){ n<-length(tr); po<-mean(tr==ca)
  pe<-mean(tr)*mean(ca)+mean(!tr)*mean(!ca); (po-pe)/(1-pe) }
boot_kappa <- function(tr,ca,B=2000){ n<-length(tr)
  v<-replicate(B,{i<-sample.int(n,n,replace=TRUE); ckappa(tr[i],ca[i])}); quantile(v,c(.025,.975),na.rm=TRUE) }
pct <- function(x,n) sprintf("%.1f%% (%.1f-%.1f)", 100*x/n, 100*wilson(x,n)[1], 100*wilson(x,n)[2])
metrics <- function(tr, ca, name){
  TP<-sum(tr&ca); FP<-sum(!tr&ca); TN<-sum(!tr&!ca); FN<-sum(tr&!ca)
  sens<-TP/(TP+FN); spec<-TN/(TN+FP); ppv<-TP/(TP+FP); npv<-TN/(TN+FN)
  acc<-(TP+TN)/length(tr); bal<-(sens+spec)/2; k<-ckappa(tr,ca); kci<-boot_kappa(tr,ca)
  data.frame(rule=name, TP,FP,TN,FN,
    sensitivity=pct(TP,TP+FN), specificity=pct(TN,TN+FP), PPV=pct(TP,TP+FP), NPV=pct(TN,TN+FN),
    accuracy=pct(TP+TN,length(tr)), bal_acc=sprintf("%.1f%%",100*bal),
    kappa=sprintf("%.3f (%.3f-%.3f)", k, kci[1], kci[2]), stringsAsFactors=FALSE)
}
rules <- list(
  "EBV RPM>=1 (PRIMARY)" = gt$ebv_rpm>=1,
  "EBV RPM>=Youden"      = gt$ebv_rpm>=thr_youden,
  "EBV RPM>=5"           = gt$ebv_rpm>=5,
  "pipeline high_conf"   = gt$ebv_high_conf=="TRUE",
  "any read >=1"         = gt$ebv_present==1)
M <- do.call(rbind, lapply(names(rules), function(n) metrics(truth, rules[[n]], n)))
M$AUC <- ""; M$AUC[1] <- sprintf("%.3f (%.3f-%.3f)", as.numeric(auc(roc1)), auc_ci[1], auc_ci[3])
write.csv(M, "table/ebv_concordance/ebv_concordance_metrics.csv", row.names=FALSE)
cat("\n================ CONCORDANCE METRICS ================\n"); print(M[,c("rule","TP","FP","TN","FN","sensitivity","specificity","kappa")], row.names=FALSE)

## ---- McNemar comparing two RULES against clinical truth (STRATIFIED, per reviewer consensus) ----
# A marginal McNemar on the two calls only tests positive-rate homogeneity, not diagnostic accuracy.
# To compare SENSITIVITY, McNemar within clinical-positives (does correctness differ on diseased);
# to compare SPECIFICITY, McNemar within clinical-negatives.
mcstrat <- function(callA, callB, lab){
  sp <- tryCatch(mcnemar.test(table(factor(callA[truth],c(F,T)),  factor(callB[truth],c(F,T))))$p.value,  error=function(e)NA)
  sn <- tryCatch(mcnemar.test(table(factor(callA[!truth],c(F,T)), factor(callB[!truth],c(F,T))))$p.value, error=function(e)NA)
  cat(sprintf("McNemar %s : sensitivity(within clin+) p=%.3g | specificity(within clin-) p=%.3g\n", lab, sp, sn))
  data.frame(comparison=lab, mcnemar_sens_p=sp, mcnemar_spec_p=sn, stringsAsFactors=FALSE)
}
cat("\n")
mcdf <- rbind(mcstrat(gt$ebv_rpm>=1, gt$ebv_high_conf=="TRUE", "RPM>=1 vs high_conf"),
              mcstrat(gt$ebv_rpm>=1, gt$ebv_present==1,        "RPM>=1 vs any-read"))
write.csv(mcdf, "table/ebv_concordance/ebv_mcnemar.csv", row.names=FALSE)
# persist Youden threshold + AUC to a machine-readable file; flag the silent 1.7 fallback
write.csv(data.frame(youden_threshold_RPM=thr_youden,
                     youden_is_hardcoded_fallback=!is.finite(as.numeric(yj$threshold[1])),
                     auc=as.numeric(auc(roc1)), auc_lo=auc_ci[1], auc_hi=auc_ci[3]),
          "table/ebv_concordance/ebv_concordance_aux.csv", row.names=FALSE)

## ---- discordants ----
disc <- gt[(truth & gt$ebv_rpm<1) | (!truth & gt$ebv_rpm>=1), c("sample","clinical_ebv","ebv_reads","ebv_rpm","ebv_high_conf")]
disc$type <- ifelse(disc$clinical_ebv=="pos","clin+ / pipeline- (microbial non-detection)",
             ifelse(disc$ebv_rpm>=5,"clin- / pipeline+ HIGH (>=5 RPM: pathology re-review)",
                    "clin- / pipeline+ low (1-5 RPM: technical QC)"))
disc <- disc[order(disc$clinical_ebv, -disc$ebv_rpm),]
write.csv(disc, "table/ebv_concordance/ebv_discordant_cases.csv", row.names=FALSE)
cat("\n================ DISCORDANT CASES ================\n"); print(disc, row.names=FALSE)
cat(sprintf("\n  pathology re-review flags (clin- & RPM>=5): %s\n",
            paste(disc$sample[grepl("re-review",disc$type)], collapse=", ")))

## ---- classify clinical-NA ----
na <- d[is.na(d$clinical_ebv), c("sample","ebv_reads","ebv_rpm","ebv_high_conf")]
na$pipeline_pred <- ifelse(na$ebv_rpm>=1,"EBV+","EBV-"); na$high_burden <- na$ebv_rpm>=5
na <- na[order(-na$ebv_rpm),]
write.csv(na, "table/ebv_concordance/ebv_NA_predictions.csv", row.names=FALSE)
cat(sprintf("\nclinical-NA (n=%d): predicted EBV+ (RPM>=1) = %d; high-burden (RPM>=5) = %d\n",
            nrow(na), sum(na$pipeline_pred=="EBV+"), sum(na$high_burden)))

## ================ FIGURES ================
# (A) ROC
rocdf <- data.frame(fpr=1-roc1$specificities, tpr=roc1$sensitivities)
pA <- ggplot(rocdf, aes(fpr,tpr)) + geom_abline(slope=1,intercept=0,linetype=2,colour="grey70") +
  geom_path(colour="#2C7FB8", linewidth=1) +
  annotate("text", .6,.2, label=sprintf("AUC = %.3f\n(DeLong 95%% CI %.3f-%.3f)",
           as.numeric(auc(roc1)),auc_ci[1],auc_ci[3]), size=4) +
  labs(title="A  ROC: EBV RPM predicts clinical EBER", x="1 - specificity", y="sensitivity") +
  coord_equal() + theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))
# (B) EBV RPM by clinical status w/ thresholds + discordants highlighted
pcps <- 0.5*1e6/median(d$ebv_reads/(pmax(d$ebv_rpm,1e-9))*1e6, na.rm=TRUE)  # not used; simple floor below
gt$grp <- ifelse(truth,"clinical EBV+","clinical EBV-")
gt$disc <- (truth & gt$ebv_rpm<1) | (!truth & gt$ebv_rpm>=1)
pB <- ggplot(gt, aes(grp, ebv_rpm+0.03)) +
  geom_hline(yintercept=c(1,5), linetype=c(2,3), colour=c("#2C7FB8","#D95F02")) +
  geom_boxplot(aes(fill=grp), outlier.shape=NA, width=.55, alpha=.5) +
  geom_jitter(aes(colour=disc), width=.15, size=1, alpha=.7) +
  scale_y_log10(labels=comma) +
  scale_fill_manual(values=c("clinical EBV+"="#E7298A","clinical EBV-"="#FBB4AE"), guide="none") +
  scale_colour_manual(values=c(`FALSE`="grey55",`TRUE`="#E7298A"), name="discordant") +
  annotate("text", 0.6, 1.3, label="RPM=1", size=3, colour="#2C7FB8", hjust=0) +
  annotate("text", 0.6, 6.5, label="RPM=5", size=3, colour="#D95F02", hjust=0) +
  labs(title="B  EBV burden by clinical status (thresholds shown)",
       subtitle="discordants in pink: clin+/low-RPM or clin-/RPM>=1", x=NULL, y="EBV RPM + 0.03 (log10)") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))
# (C) confusion matrix for primary rule
ca <- gt$ebv_rpm>=1
cm <- as.data.frame(table(clinical=factor(ifelse(truth,"EBV+","EBV-"),c("EBV+","EBV-")),
                          pipeline=factor(ifelse(ca,"EBV+","EBV-"),c("EBV+","EBV-"))))
pC <- ggplot(cm, aes(pipeline, clinical, fill=Freq)) + geom_tile(colour="white") +
  geom_text(aes(label=Freq), size=6, fontface="bold") +
  scale_fill_gradient(low="#EFF3FF", high="#2C7FB8", guide="none") +
  labs(title="C  Confusion matrix — primary rule (RPM>=1)",
       subtitle=sprintf("sens %.1f%% | spec %.1f%% | kappa %.2f | AUC %.3f",
         100*sum(truth&ca)/sum(truth), 100*sum(!truth&!ca)/sum(!truth), ckappa(truth,ca), as.numeric(auc(roc1))),
       x="pipeline (microbial EBV)", y="clinical (EBER)") +
  theme_minimal(base_size=12) + theme(plot.title=element_text(face="bold"), panel.grid=element_blank())
ggsave("analysis/ebv_concordance/EBV_concordance.pdf", (pA|pC)/pB, width=12, height=10, device=cairo_pdf)
cat("\nWROTE table/ebv_concordance/ebv_concordance_metrics.csv, ebv_discordant_cases.csv, ebv_NA_predictions.csv, figures/EBV_concordance.pdf\n")
