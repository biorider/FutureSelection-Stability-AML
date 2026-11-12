## ============================================================
## HESF-AML | 05_stability.R
## Feature-selection stability under patient-level resampling
## and varying training sample size.  (H1, H2, H3, H5)
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
source(file.path(BASE,"scripts","04_hesf_core.R"))
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
MOD  <- file.path(BASE,"results","models"); TAB <- file.path(BASE,"results","tables")
dir.create(MOD, recursive=TRUE, showWarnings=FALSE); dir.create(TAB, recursive=TRUE, showWarnings=FALSE)

X  <- readRDS(file.path(PROC,"GSE37642_GPL570_expression.rds"))
cl <- read.csv(file.path(CLIN,"GSE37642_GPL570_eligible_clinical.csv"))
X  <- X[cl$patient_id,,drop=FALSE]; tt <- cl$os_time; ee <- cl$os_event
stopifnot(identical(rownames(X), cl$patient_id))

M_OUTER  <- 20      # outer replicates per training fraction
B_INNER  <- 40      # inner bootstrap repetitions for HESF / bootstrap-enet
FRACS    <- c(0.5, 0.7, 0.9, 1.0)
THRESH   <- 0.60
METHODS  <- c("uni","lasso","enet","rsf","benet","hesf","benet30","hesf30")

CLUS <- make_hesf_cluster(18)

n <- nrow(X); p <- ncol(X); gn <- colnames(X)
all_rows <- list(); within_rows <- list(); rep_log <- list()

for (f in FRACS) {
  set.seed(SEEDS$bootstrap + round(f*1000))
  Z <- lapply(METHODS, function(m) matrix(0L, M_OUTER, p, dimnames=list(NULL,gn)))
  names(Z) <- METHODS
  wstab <- matrix(NA_real_, M_OUTER, 4, dimnames=list(NULL, c("uni","lasso","enet","rsf")))
  for (r in seq_len(M_OUTER)) {
    idx <- if (f < 1) sample.int(n, round(f*n), replace=FALSE)   # subsample
           else        sample.int(n, n, replace=TRUE)            # full-size bootstrap
    Xr <- X[idx,,drop=FALSE]; tr <- tt[idx]; er <- ee[idx]
    bb <- hesf_bootstrap(Xr, tr, er, B=B_INNER, threshold=THRESH,
                         seed=SEEDS$bootstrap + round(f*1000)*100 + r, cl=CLUS)
    sf <- single_fit_sel(Xr, tr, er, seed=SEEDS$bootstrap + r)
    if (length(bb$selected_genes))  Z$hesf   [r, bb$selected_genes]  <- 1L
    if (length(bb$boot_enet_genes)) Z$benet  [r, bb$boot_enet_genes] <- 1L
    if (length(bb$selected_topk))   Z$hesf30 [r, bb$selected_topk]   <- 1L
    if (length(bb$boot_enet_topk))  Z$benet30[r, bb$boot_enet_topk]  <- 1L
    if (!is.null(sf)) for (m in c("uni","lasso","enet","rsf"))
      if (length(sf[[m]])) Z[[m]][r, sf[[m]]] <- 1L
    for (m in c("uni","lasso","enet","rsf")) wstab[r,m] <- nogueira_stability(bb$Z[[m]])
    rep_log[[length(rep_log)+1]] <- data.frame(fraction=f, replicate=r,
      n_patients=length(unique(idx)), n_deaths=sum(er),
      inner_attempted=bb$n_attempted, inner_successful=bb$n_successful,
      n_hesf=length(bb$selected_genes), n_benet=length(bb$boot_enet_genes))
    cat(sprintf("f=%.1f rep %02d/%d  inner ok %d/%d  HESF genes %d\n",
                f, r, M_OUTER, bb$n_successful, bb$n_attempted, length(bb$selected_genes)),
        file=file.path(BASE,"logs","05_progress.txt"), append=TRUE)
  }
  saveRDS(Z, file.path(MOD, sprintf("stability_Z_frac%02d.rds", round(f*100))))
  for (m in METHODS) {
    ci <- stability_ci(Z[[m]], R=2000, seed=1)
    all_rows[[length(all_rows)+1]] <- data.frame(
      fraction=f, n_train=if (f<1) round(f*n) else n, method=m,
      nogueira=nogueira_stability(Z[[m]]), ci_lo=ci[1], ci_hi=ci[2],
      mean_jaccard=mean_jaccard(Z[[m]]),
      mean_n_selected=mean(rowSums(Z[[m]])), sd_n_selected=sd(rowSums(Z[[m]])),
      n_replicates_nonempty=sum(rowSums(Z[[m]])>0))
  }
  within_rows[[length(within_rows)+1]] <- data.frame(fraction=f,
      method=colnames(wstab), within_rep_nogueira_mean=colMeans(wstab, na.rm=TRUE),
      within_rep_nogueira_sd=apply(wstab,2,sd,na.rm=TRUE))
}
parallel::stopCluster(CLUS)

stab <- do.call(rbind, all_rows); rownames(stab) <- NULL
within<- do.call(rbind, within_rows); rownames(within) <- NULL
logs  <- do.call(rbind, rep_log)
write.csv(stab,  file.path(TAB,"TABLE2_stability_by_method_and_samplesize.csv"), row.names=FALSE)
write.csv(within,file.path(TAB,"TABLE2b_within_replicate_bootstrap_stability.csv"), row.names=FALSE)
write.csv(logs,  file.path(BASE,"results","audit","stability_replicate_log.csv"), row.names=FALSE)
print(stab)
