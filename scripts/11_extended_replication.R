## ============================================================
## 11_extended_replication.R
## Gene-level external replication across ALL validation cohorts,
## always adjusted for discovery univariate effect size, plus the
## mechanism check for the cross-platform inversion.
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
suppressPackageStartupMessages({library(survival); library(matrixStats); library(pROC)})
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
TAB  <- file.path(BASE,"results","tables"); MOD <- file.path(BASE,"results","models")

cox_score_fast <- function(X, time, event){
  o <- order(time); X <- X[o,,drop=FALSE]; time <- time[o]; event <- event[o]
  n <- nrow(X); p <- ncol(X)
  rc <- function(M) colCumsums(M[n:1,,drop=FALSE])[n:1,,drop=FALSE]
  csX <- rc(X); csX2 <- rc(X*X); S0 <- n:1
  fi <- match(time, time); ev <- which(event==1L)
  if (!length(ev)) return(list(z=rep(NA_real_,p), p=rep(NA_real_,p)))
  idx <- fi[ev]; s0 <- S0[idx]
  m1 <- csX[idx,,drop=FALSE]/s0; m2 <- csX2[idx,,drop=FALSE]/s0
  U <- colSums(X[ev,,drop=FALSE]) - colSums(m1); V <- colSums(m2 - m1*m1)
  z <- U/sqrt(V); z[!is.finite(z)] <- NA_real_
  list(z=z, p=2*pnorm(-abs(z)))
}
load_cohort <- function(xf, cf){
  X <- readRDS(file.path(PROC, xf)); cl <- read.csv(file.path(CLIN, cf))
  cl <- cl[cl$patient_id %in% rownames(X), ]
  list(X=X[cl$patient_id,,drop=FALSE], cl=cl)
}
disc <- load_cohort("GSE37642_GPL570_expression.rds","GSE37642_GPL570_eligible_clinical.csv")
VAL <- list(
  GSE37642_GPL96_97 = c("GSE37642_GPL96_97_expression.rds","GSE37642_GPL96_97_eligible_clinical.csv","platform only"),
  GSE10358          = c("GSE10358_expression.rds","GSE10358_eligible_clinical.csv","population only"),
  TCGA_LAML         = c("TCGA_LAML_expression.rds","TCGA_LAML_eligible_clinical.csv","population + platform"),
  BeatAML           = c("BeatAML_expression.rds","BeatAML_eligible_clinical.csv","population + platform"),
  GSE71014          = c("GSE71014_expression.rds","GSE71014_eligible_clinical.csv","population + platform"))

boot <- readRDS(file.path(MOD,"final_discovery_bootstrap.rds"))
st   <- boot$results
dsc  <- cox_score_fast(disc$X, disc$cl$os_time, disc$cl$os_event)
DZ   <- data.frame(gene=colnames(disc$X), absz=abs(as.numeric(dsc$z)),
                   mean_expr_disc=colMeans(disc$X), stringsAsFactors=FALSE)
DZ$dir <- sign(as.numeric(dsc$z))
DZ$consensus <- st$consensus[match(DZ$gene, st$gene)]

auc_mw <- function(s,l){ r<-rank(s); n1<-sum(l==1); n0<-sum(l==0)
                         if(!n1||!n0) NA_real_ else (sum(r[l==1])-n1*(n1+1)/2)/(n1*n0) }
rows <- list(); genes_out <- list()
for (vn in names(VAL)) {
  V <- load_cohort(VAL[[vn]][1], VAL[[vn]][2])
  g <- intersect(DZ$gene, colnames(V$X))
  sc <- cox_score_fast(V$X[, g, drop=FALSE], V$cl$os_time, V$cl$os_event)
  d <- DZ[match(g, DZ$gene), ]
  d$z_ext <- as.numeric(sc$z); d$p_ext <- as.numeric(sc$p)
  d$mean_expr_ext <- colMeans(V$X[, g, drop=FALSE])
  d <- d[is.finite(d$absz) & is.finite(d$z_ext) & is.finite(d$consensus), ]
  d$replicated <- as.integer(d$p_ext < 0.05 & sign(d$z_ext) == d$dir)
  d$cohort <- vn; genes_out[[vn]] <- d

  for (scope in c("all genes","screened genes (consensus > 0)")) {
    dd <- if (scope=="all genes") d else d[d$consensus > 0, ]
    mA <- glm(replicated ~ consensus, dd, family=binomial)
    mB <- glm(replicated ~ absz,      dd, family=binomial)
    mC <- glm(replicated ~ absz + consensus, dd, family=binomial)
    ## mechanism check: add discovery mean expression
    mD <- glm(replicated ~ absz + consensus + mean_expr_disc, dd, family=binomial)
    cC <- summary(mC)$coefficients; ciC <- confint.default(mC)
    cD <- summary(mD)$coefficients
    lrt <- function(a,b) pchisq(deviance(a)-deviance(b), df.residual(a)-df.residual(b), lower.tail=FALSE)
    rows[[length(rows)+1]] <- data.frame(
      cohort=vn, shift=VAL[[vn]][3], scope=scope, n_genes=nrow(dd),
      pct_replicated=round(100*mean(dd$replicated),2),
      cor_consensus_absz=cor(dd$consensus, dd$absz, method="spearman"),
      OR_unadj_per0.1 = exp(coef(mA)["consensus"]*0.1),
      p_unadj = summary(mA)$coefficients["consensus","Pr(>|z|)"],
      OR_adj_per0.1 = exp(cC["consensus","Estimate"]*0.1),
      adj_lo=exp(ciC["consensus",1]*0.1), adj_hi=exp(ciC["consensus",2]*0.1),
      p_adj = cC["consensus","Pr(>|z|)"], LRT_p = lrt(mB,mC),
      OR_adj_plus_expr_per0.1 = exp(cD["consensus","Estimate"]*0.1),
      p_adj_plus_expr = cD["consensus","Pr(>|z|)"],
      auc_consensus = auc_mw(dd$consensus, dd$replicated),
      auc_absz      = auc_mw(dd$absz, dd$replicated),
      stringsAsFactors=FALSE)
  }
}
R <- do.call(rbind, rows)
write.csv(R, file.path(TAB,"TABLE15_extended_replication_all_cohorts.csv"), row.names=FALSE)
saveRDS(genes_out, file.path(MOD,"gene_replication_all_cohorts.rds"))
print(R[R$scope=="all genes", c("cohort","shift","n_genes","pct_replicated","OR_unadj_per0.1",
      "OR_adj_per0.1","adj_lo","adj_hi","p_adj","OR_adj_plus_expr_per0.1","auc_consensus","auc_absz")],
      row.names=FALSE, digits=3)
