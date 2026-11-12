## ============================================================
## HESF-AML | 06_nested_cv.R
## Patient-level nested CV in the discovery cohort.
## Outer 5-fold (x3 seeds) performance; inner 4-fold ridge tuning.
## Feature selection is repeated inside every training partition.
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
source(file.path(BASE,"scripts","04_hesf_core.R"))
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
MOD  <- file.path(BASE,"results","models"); TAB  <- file.path(BASE,"results","tables")

X  <- readRDS(file.path(PROC,"GSE37642_GPL570_expression.rds"))
cl <- read.csv(file.path(CLIN,"GSE37642_GPL570_eligible_clinical.csv"))
X  <- X[cl$patient_id,,drop=FALSE]; tt <- cl$os_time; ee <- cl$os_event

OUTER_V <- 5; INNER_V <- 4; B_INNER <- 25; B_OUTER <- 100
LAMBDA  <- c(0.001, 0.01, 0.1, 1, 10, 100)   # extended: the inner optimum
                                             # hit the upper bound of the original grid
SEED_LIST <- c(20261033, 20261043, 20261053)
SEL_METHODS <- c("uni","lasso","enet","rsf","benet","hesf","benet30","hesf30")
THRESH <- 0.60

CLUS <- make_hesf_cluster(18)

## selection sets for all methods on one training partition
select_all <- function(Xtr, ttr, etr, B, seed, cl=NULL){
  bb <- hesf_bootstrap(Xtr, ttr, etr, B=B, threshold=THRESH, seed=seed, cl=cl)
  sf <- single_fit_sel(Xtr, ttr, etr, seed=seed)
  if (is.null(sf)) sf <- list(uni=character(0), lasso=character(0),
                              enet=character(0), rsf=character(0))
  list(sets=list(uni=sf$uni, lasso=sf$lasso, enet=sf$enet, rsf=sf$rsf,
                 benet=bb$boot_enet_genes, hesf=bb$selected_genes,
                 benet30=bb$boot_enet_topk, hesf30=bb$selected_topk),
       stability=bb$results, n_successful=bb$n_successful, n_attempted=bb$n_attempted)
}

## ridge Cox on a fixed gene set; training-derived standardisation
fit_ridge_model <- function(Xtr, ttr, etr, genes, lambda){
  genes <- intersect(genes, colnames(Xtr))
  if (length(genes) < 2) return(NULL)
  A  <- Xtr[, genes, drop=FALSE]
  ctr <- colMeans(A); sdv <- apply(A, 2, sd); sdv[!is.finite(sdv) | sdv==0] <- 1
  As <- scale(A, center=ctr, scale=sdv)
  fit <- try(glmnet(As, Surv(ttr, etr), family="cox", alpha=0,
                    lambda=lambda, standardize=FALSE), silent=TRUE)
  if (inherits(fit,"try-error")) return(NULL)
  lp_tr <- as.numeric(predict(fit, newx=As, type="link"))
  list(fit=fit, genes=genes, center=ctr, scale=sdv, lambda=lambda,
       baseline=breslow_baseline(ttr, etr, lp_tr), lp_train=lp_tr)
}
predict_ridge <- function(m, Xnew){
  g <- m$genes
  A <- scale(Xnew[, g, drop=FALSE], center=m$center, scale=m$scale)
  as.numeric(predict(m$fit, newx=A, type="link"))
}

res_fold <- list(); oof <- list(); fold_sets <- list(); fold_stab <- list()

for (sd_i in seq_along(SEED_LIST)) {
  seed <- SEED_LIST[sd_i]
  folds <- strat_folds(ee, OUTER_V, seed)
  for (k in seq_len(OUTER_V)) {
    tr <- which(folds != k); te <- which(folds == k)
    Xtr <- X[tr,,drop=FALSE]; Xte <- X[te,,drop=FALSE]
    ttr <- tt[tr]; etr <- ee[tr]; tte <- tt[te]; ete <- ee[te]

    ## ---- inner tuning: selection repeated inside each inner fold ----
    inner <- strat_folds(etr, INNER_V, seed + 7*k)
    inner_C <- array(NA_real_, c(INNER_V, length(SEL_METHODS), length(LAMBDA)),
                     dimnames=list(NULL, SEL_METHODS, as.character(LAMBDA)))
    for (j in seq_len(INNER_V)) {
      itr <- which(inner != j); iva <- which(inner == j)
      sel_j <- select_all(Xtr[itr,,drop=FALSE], ttr[itr], etr[itr],
                          B=B_INNER, seed=seed + 100*k + j, cl=CLUS)
      for (m in SEL_METHODS) for (li in seq_along(LAMBDA)) {
        mod <- fit_ridge_model(Xtr[itr,,drop=FALSE], ttr[itr], etr[itr],
                               sel_j$sets[[m]], LAMBDA[li])
        if (is.null(mod)) next
        inner_C[j, m, li] <- cindex(ttr[iva], etr[iva],
                                    predict_ridge(mod, Xtr[iva,,drop=FALSE]))
      }
    }
    best_lambda <- vapply(SEL_METHODS, function(m){
      mc <- colMeans(inner_C[,m,,drop=FALSE][,1,], na.rm=TRUE)
      if (all(!is.finite(mc))) LAMBDA[2] else LAMBDA[which.max(mc)] }, 0)

    ## ---- outer training: selection refit on the whole training fold ----
    sel_o <- select_all(Xtr, ttr, etr, B=B_OUTER, seed=seed + k, cl=CLUS)
    fold_sets[[length(fold_sets)+1]] <- data.frame(seed=seed, fold=k,
      method=SEL_METHODS,
      n_genes=vapply(SEL_METHODS, function(m) length(sel_o$sets[[m]]), 0L),
      best_lambda=best_lambda[SEL_METHODS],
      genes=vapply(SEL_METHODS, function(m) paste(sel_o$sets[[m]], collapse=";"), ""),
      stringsAsFactors=FALSE)
    fs <- sel_o$stability; fs$seed <- seed; fs$fold <- k
    fold_stab[[length(fold_stab)+1]] <- fs[fs$consensus > 0, ]

    for (m in SEL_METHODS) {
      mod <- fit_ridge_model(Xtr, ttr, etr, sel_o$sets[[m]], best_lambda[[m]])
      if (is.null(mod)) { res_fold[[length(res_fold)+1]] <- data.frame(seed=seed, fold=k,
        method=m, n_genes=length(sel_o$sets[[m]]), lambda=best_lambda[[m]],
        c_index=NA_real_, n_test=length(te), deaths_test=sum(ete)); next }
      lp <- predict_ridge(mod, Xte)
      res_fold[[length(res_fold)+1]] <- data.frame(seed=seed, fold=k, method=m,
        n_genes=length(mod$genes), lambda=best_lambda[[m]],
        c_index=cindex(tte, ete, lp), n_test=length(te), deaths_test=sum(ete))
      oof[[length(oof)+1]] <- data.frame(seed=seed, fold=k, method=m,
        patient_id=rownames(Xte), os_time=tte, os_event=ete, lp=lp)
    }
    ## model-class comparator: random survival forest as the predictor itself
    sfc <- screen_univariate(Xtr, ttr, etr, 0.20, 500)
    if (length(sfc) >= 2) {
      nm <- sfc; sfn <- paste0("g", seq_along(nm))
      dtr <- data.frame(.time=ttr, .event=etr, Xtr[, nm, drop=FALSE]); names(dtr) <- c(".time",".event",sfn)
      dte <- data.frame(Xte[, nm, drop=FALSE]); names(dte) <- sfn
      rf <- try(ranger(Surv(.time,.event) ~ ., data=dtr, num.trees=500,
                 mtry=max(1,floor(sqrt(length(nm)))), min.node.size=10,
                 splitrule="logrank", num.threads=18, seed=seed+k), silent=TRUE)
      if (!inherits(rf,"try-error")) {
        pr <- predict(rf, data=dte)
        risk <- rowSums(pr$chf)
        res_fold[[length(res_fold)+1]] <- data.frame(seed=seed, fold=k, method="rsf_model",
          n_genes=length(nm), lambda=NA_real_, c_index=cindex(tte, ete, risk),
          n_test=length(te), deaths_test=sum(ete))
        oof[[length(oof)+1]] <- data.frame(seed=seed, fold=k, method="rsf_model",
          patient_id=rownames(Xte), os_time=tte, os_event=ete, lp=risk)
      }
    }
    message(sprintf("seed %d fold %d done | HESF genes %d | inner ok %d/%d",
            seed, k, length(sel_o$sets$hesf), sel_o$n_successful, sel_o$n_attempted))
  }
}
parallel::stopCluster(CLUS)

RES <- do.call(rbind, res_fold); OOF <- do.call(rbind, oof)
SETS <- do.call(rbind, fold_sets); STAB <- do.call(rbind, fold_stab)
write.csv(RES,  file.path(TAB,"nested_cv_fold_results.csv"), row.names=FALSE)
write.csv(SETS, file.path(TAB,"nested_cv_selected_genes_by_fold.csv"), row.names=FALSE)
saveRDS(list(RES=RES, OOF=OOF, SETS=SETS, STAB=STAB), file.path(MOD,"nested_cv_raw.rds"))

## ---- pooled out-of-fold performance per method (per seed, then averaged) ----
summ <- do.call(rbind, lapply(split(OOF, list(OOF$method, OOF$seed), drop=TRUE), function(d){
  data.frame(method=d$method[1], seed=d$seed[1], n=nrow(d), deaths=sum(d$os_event),
             c_index=cindex(d$os_time, d$os_event, d$lp)) }))
pooled <- do.call(rbind, lapply(split(summ, summ$method), function(d)
  data.frame(method=d$method[1], mean_c_index=mean(d$c_index), sd_across_seeds=sd(d$c_index),
             min=min(d$c_index), max=max(d$c_index))))
pooled <- pooled[order(-pooled$mean_c_index), ]
write.csv(summ,   file.path(TAB,"nested_cv_pooled_by_seed.csv"), row.names=FALSE)
write.csv(pooled, file.path(TAB,"TABLE3_nested_cv_performance.csv"), row.names=FALSE)
print(pooled)
