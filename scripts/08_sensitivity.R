## ============================================================
## HESF-AML | 08_sensitivity.R
## Prespecified robustness analyses (protocol section 8).
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
source(file.path(BASE,"scripts","04_hesf_core.R"))
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

SETS <- read.csv(file.path(TAB,"nested_cv_selected_genes_by_fold.csv"))
LAM  <- median(SETS$best_lambda[SETS$method=="hesf30"])

fit_sig <- function(genes, lambda, ridge=TRUE){
  genes <- intersect(genes, colnames(X)); if (length(genes) < 2) return(NULL)
  A <- X[, genes, drop=FALSE]
  ctr <- colMeans(A); sdv <- apply(A,2,sd); sdv[!is.finite(sdv)|sdv==0] <- 1
  As <- scale(A, center=ctr, scale=sdv)
  if (ridge) { f <- glmnet(As, Surv(tt,ee), family="cox", alpha=0, lambda=lambda, standardize=FALSE)
               b <- as.matrix(coef(f))[genes,1]
  } else { d <- as.data.frame(As); f <- coxph(Surv(tt,ee) ~ ., data=d); b <- coef(f); names(b) <- genes }
  list(genes=genes, beta=b, center=ctr, scale=sdv)
}
ext_c <- function(m){
  vapply(names(val), function(vn){
    V <- val[[vn]]; g <- intersect(m$genes, colnames(V$X)); if (length(g) < 2) return(NA_real_)
    A <- scale(V$X[,g,drop=FALSE], center=m$center[g], scale=m$scale[g]); A[!is.finite(A)] <- 0
    cindex(V$cl$os_time, V$cl$os_event, as.numeric(A %*% m$beta[g])) }, 0)
}

CLUS <- make_hesf_cluster(18)
configs <- rbind(
  data.frame(label="primary (alpha=0.5, k=30, cand=500, B=100)", alpha=0.5, top_k=30, max_genes=500, B=100),
  data.frame(label="B = 200",               alpha=0.5,  top_k=30, max_genes=500,  B=200),
  data.frame(label="alpha = 0.25",          alpha=0.25, top_k=30, max_genes=500,  B=100),
  data.frame(label="alpha = 0.75",          alpha=0.75, top_k=30, max_genes=500,  B=100),
  data.frame(label="alpha = 1.00 (LASSO)",  alpha=1.00, top_k=30, max_genes=500,  B=100),
  data.frame(label="forest top-20",         alpha=0.5,  top_k=20, max_genes=500,  B=100),
  data.frame(label="forest top-50",         alpha=0.5,  top_k=50, max_genes=500,  B=100),
  data.frame(label="candidates = 200",      alpha=0.5,  top_k=30, max_genes=200,  B=100),
  data.frame(label="candidates = 1000",     alpha=0.5,  top_k=30, max_genes=1000, B=100))

W_GRID <- c(0.50, 0.60, 0.70); T_GRID <- c(0.50, 0.60, 0.70, 0.80)
size_rows <- list(); perf_rows <- list(); store <- list()
for (i in seq_len(nrow(configs))) {
  cfg <- configs[i,]
  bb <- hesf_bootstrap(X, tt, ee, B=cfg$B, alpha=cfg$alpha, top_k=cfg$top_k,
                       max_genes=cfg$max_genes, seed=SEEDS$sensitivity+i, cl=CLUS)
  store[[cfg$label]] <- bb$results
  ## threshold x weight grid is a free post-hoc recombination of the frequencies
  for (w in W_GRID) for (th in T_GRID) {
    cons <- w*bb$results$enet_frequency + (1-w)*bb$results$forest_frequency
    size_rows[[length(size_rows)+1]] <- data.frame(config=cfg$label, w_enet=w, threshold=th,
      n_genes=sum(cons >= th, na.rm=TRUE), max_consensus=round(max(cons, na.rm=TRUE),3))
  }
  sig30 <- bb$results$gene[seq_len(min(30, nrow(bb$results)))]
  for (ridge in c(TRUE, FALSE)) {
    m <- fit_sig(sig30, LAM, ridge=ridge); if (is.null(m)) next
    ec <- ext_c(m)
    perf_rows[[length(perf_rows)+1]] <- data.frame(config=cfg$label,
      survival_model=if (ridge) "ridge Cox" else "unpenalised Cox",
      n_genes=length(m$genes), n_successful_boot=bb$n_successful,
      c_GSE37642_GPL96_97=ec["GSE37642_GPL96_97"], c_GSE10358=ec["GSE10358"],
      c_TCGA_LAML=ec["TCGA_LAML"])
  }
  cat(sprintf("config %d/%d done: %s (%d/%d boots ok)\n", i, nrow(configs),
      cfg$label, bb$n_successful, bb$n_attempted),
      file=file.path(BASE,"logs","08_progress.txt"), append=TRUE)
}
parallel::stopCluster(CLUS)
saveRDS(store, file.path(MOD,"sensitivity_frequency_tables.rds"))
SZ <- do.call(rbind, size_rows); PF <- do.call(rbind, perf_rows)
write.csv(SZ, file.path(TAB,"TABLE9_sensitivity_signature_size.csv"), row.names=FALSE)
write.csv(PF, file.path(TAB,"TABLE10_sensitivity_external_cindex.csv"), row.names=FALSE)
print(SZ[SZ$w_enet==0.60,]); print(PF)
