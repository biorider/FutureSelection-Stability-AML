## ============================================================
## HESF-AML  |  03_audit.R  -- sample-level feasibility audit
## ============================================================
suppressPackageStartupMessages({library(data.table)})
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
PROC <- file.path(BASE,"data","processed"); CLIN <- file.path(BASE,"data","clinical")
AUD  <- file.path(BASE,"results","audit")

audit_survival <- function(clinical, X, cohort_name){
  req <- c("sample_id","patient_id","os_time","os_event")
  stopifnot(all(req %in% names(clinical)))
  clinical$os_time  <- as.numeric(clinical$os_time)
  clinical$os_event <- as.integer(clinical$os_event)
  stopifnot(all(is.na(clinical$os_event) | clinical$os_event %in% c(0,1)))
  dup_s <- sum(duplicated(clinical$sample_id))
  dup_p <- sum(duplicated(clinical$patient_id))
  elig  <- clinical[!is.na(clinical$patient_id) & !is.na(clinical$os_time) &
                    !is.na(clinical$os_event) & clinical$os_time > 0, , drop=FALSE]
  elig  <- elig[!duplicated(elig$patient_id), , drop=FALSE]
  elig  <- elig[elig$patient_id %in% rownames(X), , drop=FALSE]
  fwrite(elig, file.path(CLIN, paste0(cohort_name,"_eligible_clinical.csv")))
  data.frame(
    cohort = cohort_name,
    platform = clinical$platform[1],
    n_geo_samples = nrow(clinical),
    n_unique_patients = length(unique(clinical$patient_id)),
    n_duplicate_sample_ids = dup_s,
    n_patients_multiple_records = dup_p,
    n_missing_os_time = sum(is.na(clinical$os_time)),
    n_missing_os_event = sum(is.na(clinical$os_event)),
    n_os_time_nonpositive = sum(!is.na(clinical$os_time) & clinical$os_time <= 0),
    n_complete_eligible = nrow(elig),
    deaths = sum(elig$os_event==1), censored = sum(elig$os_event==0),
    median_followup_days = round(median(elig$os_time),1),
    min_followup_days = round(min(elig$os_time),1),
    max_followup_days = round(max(elig$os_time),1),
    n_genes = ncol(X),
    expression_missing = sum(is.na(X)),
    age_available = sum(!is.na(clinical$age)),
    stringsAsFactors = FALSE)
}

co <- list(
  list("GSE37642_GPL570",  "GSE37642_GPL570_clinical.csv",  "GSE37642_GPL570_expression.rds"),
  list("GSE37642_GPL96_97","GSE37642_GPL96_97_clinical.csv","GSE37642_GPL96_97_expression.rds"),
  list("GSE10358",         "GSE10358_clinical.csv",         "GSE10358_expression.rds"),
  list("TCGA_LAML",        "TCGA_LAML_clinical.csv",        "TCGA_LAML_expression.rds"))

res <- do.call(rbind, lapply(co, function(z){
  cl <- fread(file.path(CLIN, z[[2]]), data.table=FALSE)
  if (!"age" %in% names(cl)) cl$age <- NA_real_
  X  <- readRDS(file.path(PROC, z[[3]]))
  audit_survival(cl, X, z[[1]])
}))
fwrite(res, file.path(AUD,"FINAL_COHORT_FEASIBILITY.csv"))

## headline table requested in the protocol
head_tbl <- res[,c("cohort","n_unique_patients","n_complete_eligible","deaths","censored")]
names(head_tbl) <- c("Cohort","Unique patients","Complete OS","Deaths","Censored")
fwrite(head_tbl, file.path(AUD,"TABLE1_cohort_audit.csv"))

## gene-universe overlap with the discovery cohort
Xd <- readRDS(file.path(PROC,"GSE37642_GPL570_expression.rds"))
ov <- do.call(rbind, lapply(co[-1], function(z){
  X <- readRDS(file.path(PROC, z[[3]]))
  data.frame(cohort=z[[1]], n_genes=ncol(X),
             n_shared_with_discovery=length(intersect(colnames(Xd), colnames(X))),
             pct_of_discovery=round(100*length(intersect(colnames(Xd),colnames(X)))/ncol(Xd),1))
}))
fwrite(ov, file.path(AUD,"gene_universe_overlap.csv"))

## patient overlap between GSE37642 platforms (same series!)
c570 <- fread(file.path(CLIN,"GSE37642_GPL570_clinical.csv"), data.table=FALSE)
c9697<- fread(file.path(CLIN,"GSE37642_GPL96_97_clinical.csv"), data.table=FALSE)
ovp <- data.frame(n_GPL570=nrow(c570), n_GPL96_97=nrow(c9697),
                  n_shared_patient_ids=length(intersect(c570$patient_id, c9697$patient_id)))
fwrite(ovp, file.path(AUD,"GSE37642_platform_patient_overlap.csv"))

print(res[,c("cohort","n_geo_samples","n_complete_eligible","deaths","censored","median_followup_days","n_genes")])
cat("\n"); print(ov); cat("\n"); print(ovp)
