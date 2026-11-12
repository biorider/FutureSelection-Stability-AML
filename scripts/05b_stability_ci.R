## 05b_stability_ci.R -- recompute stability CIs by jackknife
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
source(file.path(BASE,"scripts","04_hesf_core.R"))
MOD <- file.path(BASE,"results","models"); TAB <- file.path(BASE,"results","tables")
st <- read.csv(file.path(TAB,"TABLE2_stability_by_method_and_samplesize.csv"))
out <- list()
for (f in unique(st$fraction)) {
  Z <- readRDS(file.path(MOD, sprintf("stability_Z_frac%02d.rds", round(f*100))))
  for (m in names(Z)) {
    ci <- stability_ci_jack(Z[[m]])
    out[[length(out)+1]] <- data.frame(fraction=f, method=m,
      nogueira=unname(ci[1]), ci_lo=unname(ci[2]), ci_hi=unname(ci[3]))
  }
}
J <- do.call(rbind, out)
st$ci_lo <- NULL; st$ci_hi <- NULL; st$nogueira <- NULL
st2 <- merge(J, st, by=c("fraction","method"))
st2$resampling <- ifelse(st2$fraction < 1, "subsample without replacement",
                                           "bootstrap (full size, with replacement)")
st2 <- st2[order(st2$fraction, -st2$nogueira), ]
write.csv(st2, file.path(TAB,"TABLE2_stability_by_method_and_samplesize.csv"), row.names=FALSE)
print(st2[,c("fraction","method","nogueira","ci_lo","ci_hi","mean_jaccard","mean_n_selected")])
