## ============================================================
## 12_revalidate_all.R
## Frozen discovery models applied to ALL five validation cohorts,
## under frozen and cohort-standardised scaling, with coverage.
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
suppressPackageStartupMessages({library(survival); library(matrixStats); library(timeROC)})
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
TAB  <- file.path(BASE,"results","tables"); MOD <- file.path(BASE,"results","models")

load_cohort <- function(xf, cf){
  X <- readRDS(file.path(PROC, xf)); cl <- read.csv(file.path(CLIN, cf))
  cl <- cl[cl$patient_id %in% rownames(X), ]
  list(X=X[cl$patient_id,,drop=FALSE], cl=cl)
}
VAL <- list(
  GSE37642_GPL96_97 = c("GSE37642_GPL96_97_expression.rds","GSE37642_GPL96_97_eligible_clinical.csv","platform only"),
  GSE10358          = c("GSE10358_expression.rds","GSE10358_eligible_clinical.csv","population only"),
  TCGA_LAML         = c("TCGA_LAML_expression.rds","TCGA_LAML_eligible_clinical.csv","population + platform"),
  BeatAML           = c("BeatAML_expression.rds","BeatAML_eligible_clinical.csv","population + platform"),
  GSE71014          = c("GSE71014_expression.rds","GSE71014_eligible_clinical.csv","population + platform"))
MODELS <- readRDS(file.path(MOD,"HESF_AML_FINAL.rds"))

cindex <- function(t,e,s) if (length(unique(s))<2) NA_real_ else
  as.numeric(concordance(Surv(t,e) ~ s, reverse=TRUE)$concordance)
breslow_baseline <- function(time,event,lp){
  o <- order(time); time<-time[o]; event<-event[o]; r<-exp(lp[o])
  S0 <- rev(cumsum(rev(r))); ti <- time[event==1]; ut <- unique(ti)
  dh <- vapply(ut, function(t) sum(ti==t)/S0[match(t,time)], 0)
  data.frame(time=ut, H0=cumsum(dh))
}
surv_at <- function(base, lp, t){
  H <- approx(base$time, base$H0, xout=t, method="constant", rule=2, f=0, ties="ordered")$y
  outer(exp(lp), H, function(a,b) exp(-b*a))
}
brier_ipcw <- function(time,event,Smat,times){
  kmC <- survfit(Surv(time, 1-event) ~ 1)
  Gf <- function(t){ s <- summary(kmC, times=pmin(t,max(kmC$time)), extend=TRUE)$surv
                     s[s<=0] <- min(s[s>0],na.rm=TRUE); s }
  vapply(seq_along(times), function(k){ t<-times[k]; S<-Smat[,k]
    mean(as.numeric(time<=t & event==1)/Gf(time)*(0-S)^2 + as.numeric(time>t)/Gf(t)*(1-S)^2) }, 0)
}
rows <- list()
for (vn in names(VAL)) {
  V <- load_cohort(VAL[[vn]][1], VAL[[vn]][2])
  for (m in names(MODELS)) {
    mod <- MODELS[[m]]; if (is.null(mod)) next
    present <- intersect(mod$genes, colnames(V$X))
    if (length(present) < 2) next
    for (mode in c("frozen","cohort_z")) {
      A <- V$X[, present, drop=FALSE]
      A <- if (mode=="frozen") scale(A, center=mod$center[present], scale=mod$scale[present]) else scale(A)
      A[!is.finite(A)] <- 0
      lp <- as.numeric(A %*% mod$beta[present])
      tt <- V$cl$os_time; ee <- V$cl$os_event
      z  <- (lp-mean(lp))/sd(lp)
      cx <- summary(coxph(Surv(tt,ee) ~ z))
      tms <- c(365,1095,1825); tms <- tms[tms < max(tt[ee==1])]
      auc <- rep(NA_real_,3)
      if (length(tms)) { tr <- try(timeROC(T=tt, delta=ee, marker=lp, cause=1, times=tms, iid=FALSE), silent=TRUE)
                         if (!inherits(tr,"try-error")) auc[seq_along(tr$AUC)] <- as.numeric(tr$AUC) }
      br <- rep(NA_real_,3)
      if (length(tms)) br[seq_along(tms)] <- brier_ipcw(tt, ee, surv_at(mod$baseline, lp, tms), tms)
      ci <- concordance(Surv(tt,ee) ~ lp, reverse=TRUE)
      rows[[length(rows)+1]] <- data.frame(
        cohort=vn, shift=VAL[[vn]][3], method=m, scaling=mode,
        n=length(tt), deaths=sum(ee),
        genes_in_signature=length(mod$genes), genes_present=length(present),
        pct_present=round(100*length(present)/length(mod$genes),1),
        c_index=as.numeric(ci$concordance),
        c_lo=as.numeric(ci$concordance)-1.96*sqrt(as.numeric(ci$var)),
        c_hi=as.numeric(ci$concordance)+1.96*sqrt(as.numeric(ci$var)),
        HR_per_SD=cx$coefficients[1,"exp(coef)"], p_value=cx$coefficients[1,"Pr(>|z|)"],
        AUC_1y=auc[1], AUC_3y=auc[2], AUC_5y=auc[3],
        Brier_1y=br[1], Brier_3y=br[2], Brier_5y=br[3], stringsAsFactors=FALSE)
    }
  }
}
P <- do.call(rbind, rows)
write.csv(P, file.path(TAB,"TABLE17_external_validation_all_cohorts.csv"), row.names=FALSE)
h <- P[P$method=="hesf30", ]
print(h[,c("cohort","shift","scaling","n","deaths","pct_present","c_index","c_lo","c_hi","HR_per_SD","p_value","Brier_3y")],
      row.names=FALSE, digits=3)
