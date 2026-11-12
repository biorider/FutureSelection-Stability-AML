## ============================================================
## 09_conditional_h4.R
## Does the HESF consensus score predict external replication
## BEYOND the discovery univariate effect size?
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
source(file.path(BASE,"scripts","04_hesf_core.R"))
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
TAB  <- file.path(BASE,"results","tables"); MOD <- file.path(BASE,"results","models")
has_proc <- requireNamespace("pROC", quietly=TRUE)

load_cohort <- function(xf, cf){
  X <- readRDS(file.path(PROC, xf)); cl <- read.csv(file.path(CLIN, cf))
  cl <- cl[cl$patient_id %in% rownames(X), ]
  list(X=X[cl$patient_id,,drop=FALSE], cl=cl)
}
disc <- load_cohort("GSE37642_GPL570_expression.rds","GSE37642_GPL570_eligible_clinical.csv")
val  <- list(
  GSE37642_GPL96_97 = load_cohort("GSE37642_GPL96_97_expression.rds","GSE37642_GPL96_97_eligible_clinical.csv"),
  GSE10358          = load_cohort("GSE10358_expression.rds","GSE10358_eligible_clinical.csv"),
  TCGA_LAML         = load_cohort("TCGA_LAML_expression.rds","TCGA_LAML_eligible_clinical.csv"))

boot <- readRDS(file.path(MOD,"final_discovery_bootstrap.rds"))
st   <- boot$results                       # gene, enet_frequency, forest_frequency, consensus
dsc  <- cox_score_fast(disc$X, disc$cl$os_time, disc$cl$os_event)
disc_z <- data.frame(gene=colnames(disc$X), z_disc=as.numeric(dsc$z), p_disc=as.numeric(dsc$p),
                     stringsAsFactors=FALSE)
disc_z$absz <- abs(disc_z$z_disc)

auc_mw <- function(score, lab){   # Mann-Whitney AUC
  r <- rank(score); n1 <- sum(lab==1); n0 <- sum(lab==0)
  if (n1==0 || n0==0) return(NA_real_)
  (sum(r[lab==1]) - n1*(n1+1)/2) / (n1*n0)
}

rows <- list(); strat <- list(); aucs <- list()
for (vn in names(val)) {
  V <- val[[vn]]
  g  <- intersect(disc_z$gene, colnames(V$X))
  sc <- cox_score_fast(V$X[, g, drop=FALSE], V$cl$os_time, V$cl$os_event)
  d <- data.frame(gene=g,
                  consensus = st$consensus[match(g, st$gene)],
                  absz      = disc_z$absz[match(g, disc_z$gene)],
                  dir_disc  = sign(disc_z$z_disc[match(g, disc_z$gene)]),
                  z_ext     = as.numeric(sc$z), p_ext = as.numeric(sc$p),
                  stringsAsFactors=FALSE)
  d <- d[is.finite(d$absz) & is.finite(d$z_ext) & is.finite(d$consensus), ]
  d$replicated <- as.integer(d$p_ext < 0.05 & sign(d$z_ext) == d$dir_disc)

  for (scope in c("all genes","screened genes (consensus > 0)")) {
    dd <- if (scope == "all genes") d else d[d$consensus > 0, ]
    if (nrow(dd) < 50) next
    mA <- glm(replicated ~ consensus,         data=dd, family=binomial)   # consensus only
    mB <- glm(replicated ~ absz,              data=dd, family=binomial)   # effect size only
    mC <- glm(replicated ~ absz + consensus,  data=dd, family=binomial)   # both
    lrt <- function(small, big){
      dev <- deviance(small) - deviance(big); df <- df.residual(small) - df.residual(big)
      c(chisq=dev, df=df, p=pchisq(dev, df, lower.tail=FALSE))
    }
    cC <- summary(mC)$coefficients
    ci <- suppressMessages(confint.default(mC))
    rows[[length(rows)+1]] <- data.frame(
      cohort=vn, scope=scope, n_genes=nrow(dd), n_replicated=sum(dd$replicated),
      cor_consensus_absz = cor(dd$consensus, dd$absz, method="spearman"),
      OR_consensus_unadj_per0.1 = exp(coef(mA)["consensus"]*0.1),
      p_consensus_unadj         = summary(mA)$coefficients["consensus","Pr(>|z|)"],
      OR_consensus_adj_per0.1   = exp(cC["consensus","Estimate"]*0.1),
      adj_lo = exp(ci["consensus",1]*0.1), adj_hi = exp(ci["consensus",2]*0.1),
      p_consensus_adj           = cC["consensus","Pr(>|z|)"],
      OR_absz_adj_per1SD        = exp(cC["absz","Estimate"]*sd(dd$absz)),
      p_absz_adj                = cC["absz","Pr(>|z|)"],
      LRT_consensus_added_chisq = lrt(mB, mC)["chisq"], LRT_consensus_added_p = lrt(mB, mC)["p"],
      LRT_absz_added_chisq      = lrt(mA, mC)["chisq"], LRT_absz_added_p      = lrt(mA, mC)["p"],
      AIC_consensus=AIC(mA), AIC_absz=AIC(mB), AIC_both=AIC(mC),
      stringsAsFactors=FALSE)

    ## discrimination for predicting replication
    aA <- auc_mw(dd$consensus, dd$replicated); aB <- auc_mw(dd$absz, dd$replicated)
    aC <- auc_mw(predict(mC, type="link"), dd$replicated)
    ar <- data.frame(cohort=vn, scope=scope, auc_consensus=aA, auc_absz=aB, auc_both=aC,
                     stringsAsFactors=FALSE)
    if (has_proc) {
      r1 <- pROC::roc(dd$replicated, dd$consensus, quiet=TRUE)
      r2 <- pROC::roc(dd$replicated, dd$absz,      quiet=TRUE)
      dl <- pROC::roc.test(r1, r2, method="delong")
      ar$delong_p_consensus_vs_absz <- dl$p.value
      ci1 <- as.numeric(pROC::ci.auc(r1)); ci2 <- as.numeric(pROC::ci.auc(r2))
      ar$auc_consensus_lo <- ci1[1]; ar$auc_consensus_hi <- ci1[3]
      ar$auc_absz_lo <- ci2[1];      ar$auc_absz_hi <- ci2[3]
    }
    aucs[[length(aucs)+1]] <- ar
  }

  ## non-parametric: replication rate by consensus tertile WITHIN quintiles of |z_disc|
  dd <- d[d$consensus > 0, ]
  dd$zq <- cut(dd$absz, breaks=quantile(dd$absz, 0:5/5), include.lowest=TRUE,
               labels=paste0("Q", 1:5))
  dd$ct <- cut(dd$consensus, breaks=unique(quantile(dd$consensus, 0:3/3)),
               include.lowest=TRUE, labels=FALSE)
  ag <- aggregate(replicated ~ zq + ct, dd, function(x) c(n=length(x), r=sum(x)))
  strat[[length(strat)+1]] <- data.frame(cohort=vn, absz_quintile=ag$zq,
    consensus_tertile=ag$ct, n_genes=ag$replicated[,"n"], n_replicated=ag$replicated[,"r"],
    pct_replicated=round(100*ag$replicated[,"r"]/ag$replicated[,"n"],1))
}
R <- do.call(rbind, rows); A <- do.call(rbind, aucs); S <- do.call(rbind, strat)
write.csv(R, file.path(TAB,"TABLE11_conditional_replication.csv"), row.names=FALSE)
write.csv(A, file.path(TAB,"TABLE12_replication_discrimination_AUC.csv"), row.names=FALSE)
write.csv(S, file.path(TAB,"TABLE13_replication_by_consensus_within_effectsize.csv"), row.names=FALSE)
cat("pROC available:", has_proc, "\n\n")
print(R[,c("cohort","scope","n_genes","cor_consensus_absz","OR_consensus_unadj_per0.1",
           "OR_consensus_adj_per0.1","adj_lo","adj_hi","p_consensus_adj","LRT_consensus_added_p")],
      row.names=FALSE, digits=3)
cat("\n"); print(A, row.names=FALSE, digits=3)
