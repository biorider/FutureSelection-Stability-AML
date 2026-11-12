## ============================================================
## HESF-AML | 07_final_external.R
## Freeze the discovery model, validate in 3 independent cohorts,
## and test whether selection stability predicts external
## reproducibility of individual genes (H4).
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
source(file.path(BASE,"scripts","04_hesf_core.R"))
suppressPackageStartupMessages(library(timeROC))
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
MOD  <- file.path(BASE,"results","models"); TAB <- file.path(BASE,"results","tables")

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
X <- disc$X; tt <- disc$cl$os_time; ee <- disc$cl$os_event

## penalty: median of the inner-tuned lambdas across outer folds (nested CV)
SETS <- read.csv(file.path(TAB,"nested_cv_selected_genes_by_fold.csv"))
lam_by_method <- tapply(SETS$best_lambda, SETS$method, median)

CLUS <- make_hesf_cluster(18)
B_FINAL <- 200
boot <- hesf_bootstrap(X, tt, ee, B=B_FINAL, threshold=0.60, seed=SEEDS$discovery, cl=CLUS)
single <- single_fit_sel(X, tt, ee, seed=SEEDS$discovery)
parallel::stopCluster(CLUS)

SETS_FINAL <- list(uni=single$uni, lasso=single$lasso, enet=single$enet, rsf=single$rsf,
                   benet=boot$boot_enet_genes, hesf=boot$selected_genes,
                   benet30=boot$boot_enet_topk, hesf30=boot$selected_topk)
write.csv(boot$results[boot$results$consensus > 0, ],
          file.path(TAB,"TABLE4_discovery_selection_frequencies.csv"), row.names=FALSE)
saveRDS(boot, file.path(MOD,"final_discovery_bootstrap.rds"))

freeze_model <- function(genes, lambda){
  genes <- intersect(genes, colnames(X)); if (length(genes) < 2) return(NULL)
  A <- X[, genes, drop=FALSE]
  ctr <- colMeans(A); sdv <- apply(A,2,sd); sdv[!is.finite(sdv)|sdv==0] <- 1
  As <- scale(A, center=ctr, scale=sdv)
  fit <- glmnet(As, Surv(tt,ee), family="cox", alpha=0, lambda=lambda, standardize=FALSE)
  lp  <- as.numeric(predict(fit, newx=As, type="link"))
  list(genes=genes, model=fit, center=ctr, scale=sdv, lambda=lambda,
       beta=as.matrix(coef(fit))[genes,1], baseline=breslow_baseline(tt,ee,lp), lp_train=lp)
}
MODELS <- lapply(names(SETS_FINAL), function(m)
  freeze_model(SETS_FINAL[[m]], unname(lam_by_method[m])))
names(MODELS) <- names(SETS_FINAL)
saveRDS(MODELS, file.path(MOD,"HESF_AML_FINAL.rds"))

predict_frozen <- function(mod, Xnew, mode=c("frozen","cohort_z")){
  mode <- match.arg(mode)
  present <- intersect(mod$genes, colnames(Xnew))
  A <- Xnew[, present, drop=FALSE]
  A <- if (mode=="frozen") scale(A, center=mod$center[present], scale=mod$scale[present])
       else                scale(A)
  A[!is.finite(A)] <- 0
  list(lp=as.numeric(A %*% mod$beta[present]), n_present=length(present),
       n_total=length(mod$genes), missing=setdiff(mod$genes, present))
}

evaluate <- function(time, event, lp, cohort, method, mode, base=NULL){
  ci <- cindex_ci(time, event, lp)
  z  <- (lp - mean(lp))/sd(lp)
  cx <- coxph(Surv(time,event) ~ z)
  sc <- summary(cx)
  ph <- tryCatch(cox.zph(cx)$table["z","p"], error=function(e) NA_real_)
  tms <- c(365,1095,1825); tms <- tms[tms < max(time[event==1])]
  auc <- rep(NA_real_, 3); auc_lo <- auc; auc_hi <- auc
  if (length(tms)) {
    tr <- try(timeROC(T=time, delta=event, marker=lp, cause=1, times=tms, iid=TRUE), silent=TRUE)
    if (!inherits(tr,"try-error")) { a <- as.numeric(tr$AUC); se <- as.numeric(tr$inference$vect_sd_1)
      auc[seq_along(a)] <- a; auc_lo[seq_along(a)] <- a-1.96*se; auc_hi[seq_along(a)] <- a+1.96*se }
  }
  br <- rep(NA_real_,3); ibs <- NA_real_
  if (!is.null(base) && length(tms)) {
    S <- surv_at(base, lp, tms); b <- brier_ipcw(time, event, S, tms)
    br[seq_along(b)] <- b
    ibs <- integrated_brier(time, event, base, lp)$ibs
  }
  data.frame(cohort=cohort, method=method, scaling=mode, n=length(time), deaths=sum(event),
    c_index=ci["est"], c_lo=ci["est"]-1.96*ci["se"], c_hi=ci["est"]+1.96*ci["se"],
    HR_per_SD=sc$coefficients[1,"exp(coef)"],
    HR_lo=sc$conf.int[1,"lower .95"], HR_hi=sc$conf.int[1,"upper .95"],
    p_value=sc$coefficients[1,"Pr(>|z|)"], ph_test_p=ph,
    AUC_1y=auc[1], AUC_1y_lo=auc_lo[1], AUC_1y_hi=auc_hi[1],
    AUC_3y=auc[2], AUC_3y_lo=auc_lo[2], AUC_3y_hi=auc_hi[2],
    AUC_5y=auc[3], AUC_5y_lo=auc_lo[3], AUC_5y_hi=auc_hi[3],
    Brier_1y=br[1], Brier_3y=br[2], Brier_5y=br[3], IBS=ibs, stringsAsFactors=FALSE)
}

rows <- list(); covg <- list(); lp_store <- list()
for (m in names(MODELS)) {
  mod <- MODELS[[m]]; if (is.null(mod)) next
  rows[[length(rows)+1]] <- evaluate(tt, ee, mod$lp_train, "GSE37642_GPL570 (apparent)",
                                     m, "frozen", mod$baseline)
  for (vn in names(val)) {
    V <- val[[vn]]
    for (mode in c("frozen","cohort_z")) {
      pr <- predict_frozen(mod, V$X, mode)
      if (pr$n_present < 2) next
      rows[[length(rows)+1]] <- evaluate(V$cl$os_time, V$cl$os_event, pr$lp, vn, m, mode, mod$baseline)
      if (mode=="frozen") {
        covg[[length(covg)+1]] <- data.frame(cohort=vn, method=m, genes_in_signature=pr$n_total,
          genes_present=pr$n_present, pct_present=round(100*pr$n_present/pr$n_total,1),
          missing_genes=paste(pr$missing, collapse=";"), stringsAsFactors=FALSE)
        lp_store[[paste(vn,m,sep="|")]] <- pr$lp
      }
    }
  }
}
## patient-level risk scores for the Kaplan-Meier / calibration figures
rs <- list()
for (vn in names(val)) {
  V <- val[[vn]]
  for (m in names(MODELS)) {
    lp <- lp_store[[paste(vn,m,sep="|")]]; if (is.null(lp)) next
    rs[[length(rs)+1]] <- data.frame(cohort=vn, method=m, patient_id=V$cl$patient_id,
      os_time=V$cl$os_time, os_event=V$cl$os_event, lp=lp,
      risk_group=ifelse(lp > median(lp), "high", "low"), stringsAsFactors=FALSE)
  }
}
for (m in names(MODELS)) {
  mod <- MODELS[[m]]; if (is.null(mod)) next
  rs[[length(rs)+1]] <- data.frame(cohort="GSE37642_GPL570 (apparent)", method=m,
    patient_id=rownames(X), os_time=tt, os_event=ee, lp=mod$lp_train,
    risk_group=ifelse(mod$lp_train > median(mod$lp_train), "high", "low"), stringsAsFactors=FALSE)
}
write.csv(do.call(rbind, rs), file.path(TAB,"risk_scores_by_cohort.csv"), row.names=FALSE)

PERF <- do.call(rbind, rows); rownames(PERF) <- NULL
write.csv(PERF, file.path(TAB,"TABLE5_external_validation_performance.csv"), row.names=FALSE)
write.csv(do.call(rbind, covg), file.path(TAB,"signature_gene_coverage_by_cohort.csv"), row.names=FALSE)

## ---- paired bootstrap: HESF vs each comparator, same patients ----
cmp <- list()
for (vn in names(val)) {
  V <- val[[vn]]
  for (ref in c("hesf","hesf30")) {
    h <- lp_store[[paste(vn,ref,sep="|")]]; if (is.null(h)) next
    for (m in setdiff(names(MODELS), ref)) {
      o <- lp_store[[paste(vn,m,sep="|")]]; if (is.null(o)) next
      b <- paired_cindex_boot(V$cl$os_time, V$cl$os_event, h, o, R=2000, seed=11)
      cmp[[length(cmp)+1]] <- data.frame(cohort=vn, comparison=paste0(ref, " - ", m),
        diff_c_index=b["diff"], lo=b["lo"], hi=b["hi"], p=b["p_two_sided"])
    }
  }
}
write.csv(do.call(rbind, cmp), file.path(TAB,"TABLE6_paired_cindex_comparisons.csv"), row.names=FALSE)

## ---- clinical adjustment (age; sex where recorded) ----
adj <- list()
for (vn in c(names(val))) {
  V <- val[[vn]]; h <- lp_store[[paste(vn,"hesf30",sep="|")]]
  if (is.null(h) || !"age" %in% names(V$cl)) next
  d <- data.frame(time=V$cl$os_time, event=V$cl$os_event,
                  z=(h-mean(h))/sd(h), age=V$cl$age)
  d <- d[complete.cases(d), ]
  if (nrow(d) < 30) next
  f1 <- coxph(Surv(time,event) ~ z, data=d); f2 <- coxph(Surv(time,event) ~ z + age, data=d)
  adj[[length(adj)+1]] <- data.frame(cohort=vn, n=nrow(d), deaths=sum(d$event),
    HR_unadj=summary(f1)$coefficients[1,"exp(coef)"], p_unadj=summary(f1)$coefficients[1,"Pr(>|z|)"],
    HR_adj_age=summary(f2)$coefficients["z","exp(coef)"],
    adj_lo=summary(f2)$conf.int["z","lower .95"], adj_hi=summary(f2)$conf.int["z","upper .95"],
    p_adj=summary(f2)$coefficients["z","Pr(>|z|)"])
}
write.csv(do.call(rbind, adj), file.path(TAB,"TABLE7_clinical_adjustment.csv"), row.names=FALSE)

## ---- H4: does discovery selection stability predict external replication? ----
st <- boot$results
disc_sc <- cox_score_fast(X, tt, ee)
disc_dir <- sign(disc_sc$z); names(disc_dir) <- colnames(X)
h4 <- list()
for (vn in names(val)) {
  V <- val[[vn]]
  g <- intersect(st$gene, colnames(V$X))
  sc <- cox_score_fast(V$X[, g, drop=FALSE], V$cl$os_time, V$cl$os_event)
  d <- data.frame(gene=g, consensus=st$consensus[match(g, st$gene)],
                  z_ext=sc$z, p_ext=sc$p, dir_disc=disc_dir[g])
  d$q_ext <- p.adjust(d$p_ext, "BH")
  d$replicated <- as.integer(d$p_ext < 0.05 & sign(d$z_ext) == d$dir_disc)
  d$bin <- cut(d$consensus, breaks=c(-Inf,0,0.2,0.4,0.6,Inf),
               labels=c("0 (never selected)","(0,0.2]","(0.2,0.4]","(0.4,0.6]",">0.6"))
  agg <- aggregate(replicated ~ bin, d, function(x) c(n=length(x), rep=sum(x)))
  out <- data.frame(cohort=vn, consensus_bin=agg$bin, n_genes=agg$replicated[,"n"],
                    n_replicated=agg$replicated[,"rep"])
  out$pct_replicated <- round(100*out$n_replicated/out$n_genes, 1)
  ## trend test: logistic regression of replication on consensus score
  lr <- glm(replicated ~ consensus, data=d, family=binomial)
  out$OR_per_0.1_consensus <- round(exp(coef(lr)["consensus"]*0.1), 3)
  out$trend_p <- signif(summary(lr)$coefficients["consensus","Pr(>|z|)"], 3)
  h4[[length(h4)+1]] <- out
  write.csv(d, file.path(TAB, paste0("H4_gene_level_replication_", vn, ".csv")), row.names=FALSE)
}
write.csv(do.call(rbind, h4), file.path(TAB,"TABLE8_H4_stability_vs_external_replication.csv"), row.names=FALSE)

cat("\n=== signature sizes ===\n"); print(sapply(SETS_FINAL, length))
cat("\n=== HESF signature (threshold 0.60) ===\n"); print(SETS_FINAL$hesf)
cat("\n=== HESF-30 signature (size-matched) ===\n"); print(SETS_FINAL$hesf30)
cat("\n=== external performance (frozen scaling) ===\n")
print(PERF[PERF$scaling=="frozen", c("cohort","method","n","deaths","c_index","c_lo","c_hi","HR_per_SD","p_value")])
