## ============================================================
## HESF-AML  |  04_hesf_core.R
## Core library: fast Cox screening, selection algorithms,
## HESF consensus, stability indices, survival metrics.
## ============================================================
suppressPackageStartupMessages({
  library(survival); library(glmnet); library(ranger)
  library(matrixStats); library(data.table)
})

SEEDS <- list(discovery=20261031, bootstrap=20261032, nested_cv=20261033,
              validation=20261034, sensitivity=20261035)

## ---------------------------------------------------------------
## Fast univariate Cox score test at beta = 0 (Breslow ties).
## Mathematically the score test reported by coxph(); computed for
## all genes at once in O(n*p) instead of p iterative fits.
## ---------------------------------------------------------------
cox_score_fast <- function(X, time, event){
  o <- order(time); X <- X[o,,drop=FALSE]; time <- time[o]; event <- event[o]
  n <- nrow(X); p <- ncol(X)
  rev_cumsum <- function(M) colCumsums(M[n:1,,drop=FALSE])[n:1,,drop=FALSE]
  csX  <- rev_cumsum(X)
  csX2 <- rev_cumsum(X*X)
  S0   <- n:1
  first_idx <- match(time, time)          # risk set start for tied times
  ev <- which(event == 1L)
  if (!length(ev)) return(list(z=rep(NA_real_,p), p=rep(NA_real_,p)))
  idx <- first_idx[ev]
  s0  <- S0[idx]
  m1  <- csX[idx,,drop=FALSE]  / s0
  m2  <- csX2[idx,,drop=FALSE] / s0
  U <- colSums(X[ev,,drop=FALSE]) - colSums(m1)
  V <- colSums(m2 - m1*m1)
  z <- U/sqrt(V); z[!is.finite(z)] <- NA_real_
  list(z=z, p=2*pnorm(-abs(z)), U=U, V=V)
}

## Component A -- univariate survival screening
screen_univariate <- function(X, time, event, p_thresh=0.20, max_genes=500){
  st <- cox_score_fast(X, time, event)
  p  <- st$p; p[!is.finite(p)] <- 1
  cand <- which(p < p_thresh)
  if (!length(cand)) return(character(0))
  cand <- cand[order(p[cand])]
  colnames(X)[cand[seq_len(min(max_genes, length(cand)))]]
}

## Component B -- penalised Cox selection (alpha = 1 LASSO, 0.5 elastic net)
fit_penalised_sel <- function(Xc, time, event, alpha){
  if (ncol(Xc) < 2) return(character(0))
  y <- Surv(time, event)
  fit <- try(cv.glmnet(Xc, y, family="cox", alpha=alpha,
                       nfolds=5, type.measure="C", standardize=TRUE), silent=TRUE)
  if (inherits(fit,"try-error")) return(NULL)
  b <- as.matrix(coef(fit, s="lambda.1se"))
  g <- rownames(b)[b[,1] != 0]
  setdiff(g, "(Intercept)")
}

## Component C -- random survival forest permutation importance
fit_rsf_sel <- function(Xc, time, event, top_k=30, seed=1L, num.threads=1){
  if (ncol(Xc) < 2) return(character(0))
  nm <- colnames(Xc); sf <- paste0("g", seq_along(nm))
  d <- data.frame(.time=time, .event=event, Xc, check.names=FALSE)
  names(d) <- c(".time",".event", sf)
  fit <- try(ranger(Surv(.time,.event) ~ ., data=d, num.trees=500,
                    mtry=max(1, floor(sqrt(ncol(Xc)))), min.node.size=10,
                    importance="permutation", splitrule="logrank",
                    num.threads=num.threads, seed=seed), silent=TRUE)
  if (inherits(fit,"try-error")) return(NULL)
  imp <- sort(fit$variable.importance, decreasing=TRUE)
  nm[match(names(imp)[seq_len(min(top_k, length(imp)))], sf)]
}

## ---------------------------------------------------------------
## One bootstrap repetition: all selection algorithms on a common
## candidate universe -> paired per-repetition selections
## ---------------------------------------------------------------
hesf_rep <- function(b, ids, X, time, event, prm){
  eb <- event[ids]
  if (sum(eb) < prm$min_deaths || length(unique(ids)) < prm$min_unique)
    return(list(ok=FALSE, why="insufficient_events_or_patients"))
  xb <- X[ids,,drop=FALSE]; tb <- time[ids]
  v <- colVars(xb); keep <- is.finite(v) & v > 0
  if (sum(keep) < 2) return(list(ok=FALSE, why="no_variable_genes"))
  xb <- xb[, keep, drop=FALSE]
  cand <- screen_univariate(xb, tb, eb, prm$p_thresh, prm$max_genes)
  if (length(cand) < 2) return(list(ok=FALSE, why="screen_empty"))
  xc <- xb[, cand, drop=FALSE]
  sel <- list(uni = head(cand, prm$top_k))
  if ("lasso" %in% prm$methods) sel$lasso <- fit_penalised_sel(xc, tb, eb, 1.0)
  if ("enet"  %in% prm$methods) sel$enet  <- fit_penalised_sel(xc, tb, eb, prm$alpha)
  if ("rsf"   %in% prm$methods) sel$rsf   <- fit_rsf_sel(xc, tb, eb, prm$top_k,
                                                seed=prm$seed+b, num.threads=1)
  if (any(vapply(sel, is.null, TRUE))) return(list(ok=FALSE, why="model_fit_failed"))
  list(ok=TRUE, sel=sel)
}

hesf_bootstrap <- function(X, time, event, B=100, alpha=0.5, p_thresh=0.20,
                           max_genes=500, top_k=30, w_enet=0.60, threshold=0.60,
                           min_deaths=10, min_unique=30, seed=20261032,
                           methods=c("uni","lasso","enet","rsf"), cl=NULL){
  force(X); force(time); force(event)   # force before cluster serialisation
  prm <- list(alpha=alpha, p_thresh=p_thresh, max_genes=max_genes, top_k=top_k,
              min_deaths=min_deaths, min_unique=min_unique, seed=seed, methods=methods)
  n <- nrow(X); p <- ncol(X); gnames <- colnames(X)
  set.seed(seed)
  IDS <- lapply(seq_len(B), function(b) sample.int(n, n, replace=TRUE))
  runner <- function(b) hesf_rep(b, IDS[[b]], X, time, event, prm)
  out <- if (is.null(cl)) lapply(seq_len(B), runner)
         else parallel::parLapply(cl, seq_len(B), runner)
  ok <- vapply(out, function(z) isTRUE(z$ok), TRUE)
  Z <- lapply(methods, function(m){
    M <- matrix(0L, nrow=sum(ok), ncol=p, dimnames=list(NULL, gnames))
    k <- 0L
    for (i in which(ok)) { k <- k+1L; g <- out[[i]]$sel[[m]]; if (length(g)) M[k, g] <- 1L }
    M })
  names(Z) <- methods
  nsel <- if (sum(ok)) sapply(Z, rowSums) else NULL
  fq <- function(nm) if (sum(ok) && !is.null(Z[[nm]])) colMeans(Z[[nm]]) else rep(NA_real_, p)
  enet_freq <- fq("enet"); rsf_freq <- fq("rsf"); uni_freq <- fq("uni")
  consensus <- w_enet*enet_freq + (1-w_enet)*rsf_freq
  res <- data.frame(gene=gnames, uni_frequency=uni_freq, enet_frequency=enet_freq,
                    forest_frequency=rsf_freq, consensus=consensus, stringsAsFactors=FALSE)
  ord_key <- if (all(is.na(consensus))) -uni_freq else -consensus
  res <- res[order(ord_key, res$gene), ]
  list(results=res, Z=Z, n_attempted=B, n_successful=sum(ok),
       failures=vapply(out[!ok], function(z) z$why, ""), n_selected=nsel,
       selected_genes  = res$gene[res$consensus      >= threshold],
       boot_enet_genes = res$gene[res$enet_frequency >= threshold],
       ## size-matched variants (same k as the RSF / univariate comparators)
       selected_topk   = res$gene[seq_len(min(top_k, nrow(res)))],
       boot_enet_topk  = res$gene[order(-res$enet_frequency)][seq_len(min(top_k, nrow(res)))],
       ## bootstrap-univariate comparator: completes the factorial of
       ## {single fit, bootstrap} x {univariate, penalised, forest, hybrid}
       boot_uni_topk   = res$gene[order(-res$uni_frequency)][seq_len(min(top_k, nrow(res)))],
       boot_uni_genes  = res$gene[!is.na(res$uni_frequency) & res$uni_frequency >= threshold],
       params=c(prm, list(w_enet=w_enet, threshold=threshold)))
}

## cluster helper -------------------------------------------------
make_hesf_cluster <- function(ncores=NULL, rng_seed=20261031){
  if (is.null(ncores)) ncores <- max(1, min(12, parallel::detectCores()-2))
  cl <- parallel::makePSOCKcluster(ncores)
  ## cv.glmnet draws its own CV folds from the worker's RNG; without
  ## L'Ecuyer streams the selected gene set is not bit-reproducible.
  parallel::clusterSetRNGStream(cl, rng_seed)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({library(survival); library(glmnet)
      library(ranger); library(matrixStats)}); NULL })
  parallel::clusterExport(cl, c("cox_score_fast","screen_univariate",
      "fit_penalised_sel","fit_rsf_sel","hesf_rep"), envir=globalenv())
  cl
}

## single (non-bootstrap) fit of each comparator on a given sample
single_fit_sel <- function(X, time, event, alpha=0.5, p_thresh=0.20, max_genes=500,
                           top_k=30, seed=1L, num.threads=1){
  v <- colVars(X); X <- X[, is.finite(v) & v>0, drop=FALSE]
  cand <- screen_univariate(X, time, event, p_thresh, max_genes)
  if (length(cand) < 2) return(NULL)
  xc <- X[, cand, drop=FALSE]
  list(uni   = head(cand, top_k),
       lasso = fit_penalised_sel(xc, time, event, 1.0),
       enet  = fit_penalised_sel(xc, time, event, alpha),
       rsf   = fit_rsf_sel(xc, time, event, top_k, seed=seed, num.threads=num.threads))
}

## ---------------------------------------------------------------
## Stability indices
## ---------------------------------------------------------------
nogueira_stability <- function(Z){
  Z <- as.matrix(Z); M <- nrow(Z); p <- ncol(Z)
  if (M < 2 || p < 2) return(NA_real_)
  k_bar <- mean(rowSums(Z))
  den <- (k_bar/p)*(1 - k_bar/p)
  if (!is.finite(den) || den == 0) return(NA_real_)
  1 - mean(colVars(Z*1.0))/den
}
mean_jaccard <- function(Z){
  Z <- as.matrix(Z); M <- nrow(Z)
  if (M < 2) return(NA_real_)
  pr <- combn(M, 2)
  sc <- apply(pr, 2, function(ix){
    a <- which(Z[ix[1],]==1); b <- which(Z[ix[2],]==1)
    d <- length(union(a,b)); if (d==0) NA_real_ else length(intersect(a,b))/d })
  mean(sc, na.rm=TRUE)
}
## Jackknife (leave-one-replicate-out) CI for the Nogueira index.
## Preferred over a row-bootstrap CI: resampling replicate rows with
## replacement duplicates rows, shrinks the per-feature variance and
## biases the index upward, so the bootstrap interval need not contain
## the point estimate.
stability_ci_jack <- function(Z){
  Z <- as.matrix(Z); M <- nrow(Z)
  if (M < 3) return(c(NA_real_, NA_real_, NA_real_))
  th  <- nogueira_stability(Z)
  thi <- vapply(seq_len(M), function(i) nogueira_stability(Z[-i,,drop=FALSE]), 0)
  se  <- sqrt((M-1)/M * sum((thi - mean(thi))^2))
  c(est=th, lo=th-1.96*se, hi=th+1.96*se)
}
stability_ci <- function(Z, R=1000, seed=1){
  set.seed(seed); M <- nrow(Z)
  if (M < 3) return(c(NA_real_, NA_real_))
  v <- replicate(R, nogueira_stability(Z[sample.int(M, M, replace=TRUE),,drop=FALSE]))
  unname(quantile(v, c(.025,.975), na.rm=TRUE))
}

## ---------------------------------------------------------------
## Survival metrics
## ---------------------------------------------------------------
cindex <- function(time, event, score){
  if (length(unique(score)) < 2) return(NA_real_)
  as.numeric(concordance(Surv(time, event) ~ score, reverse=TRUE)$concordance)
}
cindex_ci <- function(time, event, score){
  f <- concordance(Surv(time,event) ~ score, reverse=TRUE)
  c(est=as.numeric(f$concordance), se=sqrt(as.numeric(f$var)))
}
## Breslow baseline cumulative hazard from training lp
breslow_baseline <- function(time, event, lp){
  o <- order(time); time <- time[o]; event <- event[o]; r <- exp(lp[o])
  n <- length(time); S0 <- rev(cumsum(rev(r)))
  ti <- time[event==1]; idx <- match(ti, time)
  ut <- unique(ti)
  dh <- vapply(ut, function(t){ sum(ti == t) / S0[match(t, time)] }, 0)
  data.frame(time=ut, H0=cumsum(dh))
}
surv_at <- function(base, lp, t){
  H <- approx(base$time, base$H0, xout=t, method="constant", rule=2, f=0, ties="ordered")$y
  H <- ifelse(t < min(base$time), 0, H)
  outer(exp(lp), H, function(a,b) exp(-b*a))     # n x length(t)
}
## IPCW Brier score (KM censoring distribution estimated in the evaluation set)
brier_ipcw <- function(time, event, Smat, times){
  n <- length(time)
  kmC <- survfit(Surv(time, 1-event) ~ 1)
  Gf  <- function(t) { s <- summary(kmC, times=pmin(t, max(kmC$time)), extend=TRUE)$surv
                       s[s <= 0] <- min(s[s>0], na.rm=TRUE); s }
  vapply(seq_along(times), function(k){
    t <- times[k]; S <- Smat[,k]
    Gt <- Gf(t); Gi <- Gf(time)
    w1 <- as.numeric(time <= t & event==1) / Gi
    w2 <- as.numeric(time >  t) / Gt
    mean(w1*(0-S)^2 + w2*(1-S)^2)
  }, 0)
}
integrated_brier <- function(time, event, base, lp, tmax=NULL, grid=100){
  if (is.null(tmax)) tmax <- quantile(time[event==1], 0.9)
  ts <- seq(min(time[time>0])*1.01, tmax, length.out=grid)
  S  <- surv_at(base, lp, ts)
  bs <- brier_ipcw(time, event, S, ts)
  list(times=ts, brier=bs, ibs=sum(diff(ts)*(head(bs,-1)+tail(bs,-1))/2)/(max(ts)-min(ts)))
}
## paired bootstrap comparison of two C-indices (patients resampled)
paired_cindex_boot <- function(time, event, s1, s2, R=2000, seed=7){
  set.seed(seed); n <- length(time)
  d <- replicate(R, { i <- sample.int(n,n,replace=TRUE)
    if (sum(event[i]) < 5) return(NA_real_)
    cindex(time[i],event[i],s1[i]) - cindex(time[i],event[i],s2[i]) })
  d <- d[is.finite(d)]
  obs <- cindex(time,event,s1) - cindex(time,event,s2)
  c(diff=obs, lo=unname(quantile(d,.025)), hi=unname(quantile(d,.975)),
    p_two_sided=2*min(mean(d<=0), mean(d>=0)))
}

## stratified k-fold assignment on the event indicator
strat_folds <- function(event, v, seed){
  set.seed(seed); f <- integer(length(event))
  for (lev in unique(event)) { i <- which(event==lev)
    f[i] <- sample(rep_len(seq_len(v), length(i))) }
  f
}
