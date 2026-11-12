## ============================================================
## 10_null_and_swap.R
##  (a) bootstrap-univariate comparator for the discovery design
##  (b) permutation-null floor for the Nogueira index
##  (c) swapped-discovery stability curve (GPL96/97, n = 417)
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
source(file.path(BASE,"scripts","04_hesf_core.R"))
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
MOD  <- file.path(BASE,"results","models"); TAB <- file.path(BASE,"results","tables")
LOG  <- file.path(BASE,"logs","10_progress.txt")
say  <- function(...) cat(sprintf(...), file=LOG, append=TRUE)

load_cohort <- function(xf, cf){
  X <- readRDS(file.path(PROC, xf)); cl <- read.csv(file.path(CLIN, cf))
  cl <- cl[cl$patient_id %in% rownames(X), ]
  list(X=X[cl$patient_id,,drop=FALSE], cl=cl)
}
DISC <- load_cohort("GSE37642_GPL570_expression.rds","GSE37642_GPL570_eligible_clinical.csv")
SWAP <- load_cohort("GSE37642_GPL96_97_expression.rds","GSE37642_GPL96_97_eligible_clinical.csv")

FRACS  <- c(0.5, 0.7, 0.9, 1.0)
THRESH <- 0.60
CLUS   <- make_hesf_cluster(18)

## ---- generic stability experiment -----------------------------------
run_experiment <- function(X, time, event, label, methods, M, B, permute=FALSE, seed0){
  n <- nrow(X); p <- ncol(X); gn <- colnames(X)
  out_rows <- list(); Zstore <- list()
  for (f in FRACS) {
    set.seed(seed0 + round(f*1000))
    if (permute) {                       # break the gene-outcome link once per design
      o <- sample.int(n); time_u <- time[o]; event_u <- event[o]
    } else { time_u <- time; event_u <- event }
    Z <- lapply(methods, function(m) matrix(0L, M, p, dimnames=list(NULL, gn)))
    names(Z) <- methods
    for (r in seq_len(M)) {
      idx <- if (f < 1) sample.int(n, round(f*n), replace=FALSE) else sample.int(n, n, replace=TRUE)
      Xr <- X[idx,,drop=FALSE]; tr <- time_u[idx]; er <- event_u[idx]
      bb <- hesf_bootstrap(Xr, tr, er, B=B, threshold=THRESH,
                           methods=if ("enet" %in% methods || "hesf30" %in% methods)
                                     c("uni","lasso","enet","rsf") else c("uni"),
                           seed=seed0 + round(f*1000)*100 + r, cl=CLUS)
      ## NOTE: must be <<- . With `<-` the subset assignment creates a LOCAL
      ## copy of Z inside the closure and every selection matrix stays empty.
      put <- function(nm, g) if (!is.null(g) && length(g)) Z[[nm]][r, g] <<- 1L
      if ("buni30" %in% methods) put("buni30", bb$boot_uni_topk)
      if ("hesf30" %in% methods) put("hesf30", bb$selected_topk)
      if ("benet30"%in% methods) put("benet30", bb$boot_enet_topk)
      if ("hesf"   %in% methods) put("hesf",   bb$selected_genes)
      if ("benet"  %in% methods) put("benet",  bb$boot_enet_genes)
      if (any(c("uni","lasso","enet","rsf") %in% methods)) {
        sf <- single_fit_sel(Xr, tr, er, seed=seed0 + r)
        if (!is.null(sf)) for (m in intersect(c("uni","lasso","enet","rsf"), methods)) put(m, sf[[m]])
      }
      say("%s f=%.1f rep %02d/%d (n=%d, deaths=%d) inner %d/%d\n",
          label, f, r, M, length(unique(idx)), sum(er), bb$n_successful, bb$n_attempted)
    }
    filled <- vapply(Z, function(M) sum(M), 0)
    if (any(filled == 0))
      stop("empty selection matrix for: ", paste(names(filled)[filled==0], collapse=", "),
           "  (design ", label, ", fraction ", f, ")")
    Zstore[[as.character(f)]] <- Z
    for (m in methods) {
      ci <- stability_ci_jack(Z[[m]])
      out_rows[[length(out_rows)+1]] <- data.frame(
        design=label, fraction=f, n_train=if (f<1) round(f*n) else n,
        n_deaths_train=round(if (f<1) f*sum(event_u) else sum(event_u)),
        method=m, nogueira=unname(ci[1]), ci_lo=unname(ci[2]), ci_hi=unname(ci[3]),
        mean_jaccard=mean_jaccard(Z[[m]]), mean_n_selected=mean(rowSums(Z[[m]])),
        n_replicates=M, inner_B=B, permuted=permute, stringsAsFactors=FALSE)
    }
  }
  saveRDS(Zstore, file.path(MOD, paste0("stabilityZ_", label, ".rds")))
  do.call(rbind, out_rows)
}

ALL  <- c("uni","lasso","enet","rsf","benet","benet30","hesf","hesf30","buni30")
res <- list()

## (a) bootstrap-univariate comparator on the discovery design (univariate only -> fast)
res$a <- run_experiment(DISC$X, DISC$cl$os_time, DISC$cl$os_event,
                        "discovery_n136_bootuni", methods=c("buni30"),
                        M=20, B=40, permute=FALSE, seed0=SEEDS$bootstrap)
write.csv(do.call(rbind,res), file.path(TAB,"TABLE14_stability_extended.csv"), row.names=FALSE)

## (b) permutation null floors
res$b <- run_experiment(DISC$X, DISC$cl$os_time, DISC$cl$os_event,
                        "null_discovery_n136", methods=ALL, M=10, B=25,
                        permute=TRUE, seed0=SEEDS$sensitivity)
write.csv(do.call(rbind,res), file.path(TAB,"TABLE14_stability_extended.csv"), row.names=FALSE)

res$c <- run_experiment(SWAP$X, SWAP$cl$os_time, SWAP$cl$os_event,
                        "null_swapped_n417", methods=ALL, M=10, B=25,
                        permute=TRUE, seed0=SEEDS$sensitivity+1)
write.csv(do.call(rbind,res), file.path(TAB,"TABLE14_stability_extended.csv"), row.names=FALSE)

## (c) swapped-discovery stability curve
res$d <- run_experiment(SWAP$X, SWAP$cl$os_time, SWAP$cl$os_event,
                        "swapped_n417", methods=ALL, M=12, B=40,
                        permute=FALSE, seed0=SEEDS$bootstrap+7)
parallel::stopCluster(CLUS)

OUT <- do.call(rbind, res); rownames(OUT) <- NULL
write.csv(OUT, file.path(TAB,"TABLE14_stability_extended.csv"), row.names=FALSE)
print(OUT[OUT$design=="swapped_n417", c("fraction","n_train","n_deaths_train","method","nogueira","ci_lo","ci_hi")],
      row.names=FALSE, digits=3)
