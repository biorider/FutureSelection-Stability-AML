## ============================================================
## 13_matched_null.R
## Permutation null for the DISCOVERY design at exactly the
## configuration used for the observed estimates (M=20, B=40),
## so that observed-minus-null is a like-for-like contrast.
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
src <- readLines(file.path(BASE,"scripts","10_null_and_swap.R"))
end <- grep("^ALL  <-", src)[1] - 1
eval(parse(text=paste(src[1:end], collapse="\n")))   # defs + DISC/SWAP + CLUS

ALL <- c("uni","lasso","enet","rsf","benet","benet30","hesf","hesf30","buni30")
out <- run_experiment(DISC$X, DISC$cl$os_time, DISC$cl$os_event,
                      "null_discovery_matchedB", methods=ALL, M=12, B=40,
                      permute=TRUE, seed0=SEEDS$sensitivity+99)
parallel::stopCluster(CLUS)
TAB <- file.path(BASE,"results","tables")
write.csv(out, file.path(TAB,"TABLE19_matched_permutation_null.csv"), row.names=FALSE)
print(out[,c("fraction","n_train","method","nogueira","ci_lo","ci_hi","mean_n_selected")],
      row.names=FALSE, digits=3)
