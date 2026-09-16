#!/usr/bin/env Rscript
# ROC (AUC) of microbial EBV RPM vs metadata clinical EBV, with the RPM=1.7 (Youden)
# operating point marked. AUC is threshold-independent; 1.7 is one point on the curve.
suppressMessages({library(ggplot2); library(pROC)})
d  <- read.csv("data/ebv_concordance_input.csv", stringsAsFactors=FALSE)
gt <- d[d$clinical_ebv %in% c("pos","neg"), ]
truth <- gt$clinical_ebv=="pos"; score <- gt$ebv_rpm
roc1 <- roc(truth, score, quiet=TRUE, direction="<")
aci  <- ci.auc(roc1, method="delong")
rocdf <- data.frame(fpr=1-roc1$specificities, tpr=roc1$sensitivities)

# operating point at any RPM cutoff
op <- function(t){ ca<-score>=t
  data.frame(thr=t, fpr=1-sum(!truth&!ca)/sum(!truth), tpr=sum(truth&ca)/sum(truth),
             sens=sum(truth&ca)/sum(truth), spec=sum(!truth&!ca)/sum(!truth),
             TP=sum(truth&ca),FP=sum(!truth&ca),TN=sum(!truth&!ca),FN=sum(truth&!ca)) }
pts <- do.call(rbind, lapply(c(1,1.7,5), op)); pts$lab <- sprintf("RPM≥%.1f", pts$thr)
# manual label offsets so the near-coincident 1.0/1.7 points don't collide
pts$lx <- pts$fpr + c(0.03, -0.015, 0.03); pts$ly <- pts$tpr + c(0.045, -0.05, 0.0)
pts$hj <- c(0, 1, 0)
y17 <- op(1.7)

p <- ggplot(rocdf, aes(fpr,tpr)) +
  geom_abline(slope=1,intercept=0,linetype=2,colour="grey70") +
  geom_path(colour="#2C7FB8", linewidth=1) +
  geom_point(data=pts, aes(fpr,tpr), colour="grey40", size=2) +
  geom_point(data=y17, aes(fpr,tpr), colour="#E7298A", size=4) +
  annotate("segment", x=y17$fpr, xend=y17$fpr, y=0, yend=y17$tpr, linetype=3, colour="#E7298A") +
  annotate("segment", x=0, xend=y17$fpr, y=y17$tpr, yend=y17$tpr, linetype=3, colour="#E7298A") +
  geom_text(data=pts, aes(lx,ly,label=lab,hjust=hj), size=3.3,
            colour=ifelse(pts$thr==1.7,"#E7298A","grey40"), fontface=ifelse(pts$thr==1.7,"bold","plain")) +
  annotate("text", .62,.30, hjust=0, size=4.2,
    label=sprintf("AUC = %.3f\n(DeLong 95%% CI %.3f-%.3f)", as.numeric(auc(roc1)), aci[1], aci[3])) +
  annotate("text", .62,.14, hjust=0, size=3.6, colour="#E7298A",
    label=sprintf("cutoff RPM≥1.7 (Youden):\nsens %.1f%%  spec %.1f%%\nTP%d FP%d TN%d FN%d",
                  100*y17$sens,100*y17$spec,y17$TP,y17$FP,y17$TN,y17$FN)) +
  scale_x_continuous(expand=expansion(mult=c(0,.02))) + scale_y_continuous(expand=expansion(mult=c(0,.02))) +
  coord_equal() +
  labs(title="ROC: microbial EBV RPM vs clinical (metadata) EBV status",
       subtitle="cHL n=167 (75 EBER+ / 92 EBER−); operating point at RPM≥1.7 highlighted",
       x="1 - specificity (false-positive rate)", y="sensitivity (true-positive rate)") +
  theme_classic(base_size=12) + theme(plot.title=element_text(face="bold"))
ggsave("analysis/ebv_concordance/EBV_ROC_youden1.7.pdf", p, width=7.5, height=7, device=cairo_pdf)
cat(sprintf("AUC=%.3f  | RPM>=1.7: sens=%.1f%% spec=%.1f%% (TP%d FP%d TN%d FN%d)\n",
    as.numeric(auc(roc1)),100*y17$sens,100*y17$spec,y17$TP,y17$FP,y17$TN,y17$FN))
cat("WROTE analysis/ebv_concordance/EBV_ROC_youden1.7.pdf\n")
