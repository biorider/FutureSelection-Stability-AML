## ============================================================
## HESF-AML  |  02_build_tcga.R   (TCGA-LAML, RNA-seq)
## ============================================================
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
RAW  <- file.path(BASE,"data","raw"); PROC <- file.path(BASE,"data","processed")
CLIN <- file.path(BASE,"data","clinical"); AUD <- file.path(BASE,"results","audit")

ex <- file.path(RAW,"TCGA_STAR")
if (!dir.exists(ex) || length(list.files(ex, recursive=TRUE, pattern="tsv$")) < 100) {
  dir.create(ex, showWarnings=FALSE)
  untar(file.path(RAW,"TCGA_LAML_STAR_counts.tar.gz"), exdir=ex)
}
files <- fromJSON(file.path(RAW,"TCGA_LAML_files.json"), simplifyDataFrame=FALSE)
fmap <- data.frame(
  file_id   = vapply(files, function(x) x$file_id, ""),
  file_name = vapply(files, function(x) x$file_name, ""),
  case_id   = vapply(files, function(x) x$cases[[1]]$submitter_id, ""),
  sample_id = vapply(files, function(x) x$cases[[1]]$samples[[1]]$submitter_id, ""),
  stringsAsFactors = FALSE)

read_star <- function(fid, fname){
  p <- file.path(ex, fid, fname)
  d <- fread(p, sep="\t", skip=1, header=TRUE, showProgress=FALSE, data.table=FALSE)
  d <- d[!startsWith(d$gene_id, "N_"), ]
  d[d$gene_type == "protein_coding", c("gene_name","tpm_unstranded")]
}
first <- read_star(fmap$file_id[1], fmap$file_name[1])
genes <- first$gene_name
M <- matrix(NA_real_, nrow=nrow(first), ncol=nrow(fmap),
            dimnames=list(genes, fmap$case_id))
M[,1] <- first$tpm_unstranded
for (i in 2:nrow(fmap)) {
  d <- read_star(fmap$file_id[i], fmap$file_name[i])
  stopifnot(identical(d$gene_name, genes))
  M[,i] <- d$tpm_unstranded
}
## duplicate gene symbols -> keep highest-mean row (same rule as microarrays)
mu  <- rowMeans(M, na.rm=TRUE)
ord <- order(rownames(M), -mu); M <- M[ord, ]
M   <- M[!duplicated(rownames(M)), ]
Xt  <- t(log2(M + 1))                       # patients x genes, log2(TPM+1)

## ---------- clinical ----------
cs <- fromJSON(file.path(CLIN,"TCGA_LAML_cases_raw.json"), simplifyDataFrame=FALSE)
pick <- function(x, path, default=NA) { for (p in path) { if (is.null(x)) return(default); x <- x[[p]] }
                                        if (is.null(x)) default else x }
cl <- data.frame(
  patient_id = vapply(cs, function(x) pick(x,"submitter_id",NA_character_), ""),
  vital      = vapply(cs, function(x) as.character(pick(x,c("demographic","vital_status"),NA)), ""),
  d_death    = vapply(cs, function(x) as.numeric(pick(x,c("demographic","days_to_death"),NA)), 0),
  d_fup      = vapply(cs, function(x) { v <- pick(x,"diagnoses",NULL)
                      if (is.null(v)) NA_real_ else as.numeric(pick(v[[1]],"days_to_last_follow_up",NA)) }, 0),
  sex        = vapply(cs, function(x) as.character(pick(x,c("demographic","gender"),NA)), ""),
  age_days   = vapply(cs, function(x) { v <- pick(x,"diagnoses",NULL)
                      if (is.null(v)) NA_real_ else as.numeric(pick(v[[1]],"age_at_diagnosis",NA)) }, 0),
  stringsAsFactors = FALSE)
cl$os_event <- ifelse(tolower(cl$vital)=="dead",1L, ifelse(tolower(cl$vital)=="alive",0L,NA_integer_))
cl$os_time  <- ifelse(cl$os_event==1L, cl$d_death, cl$d_fup)
cl$age      <- cl$age_days/365.25
cl$sample_id <- cl$patient_id; cl$cohort <- "TCGA-LAML"; cl$platform <- "RNA-seq (STAR/TPM)"

keep <- intersect(rownames(Xt), cl$patient_id)
Xt   <- Xt[keep, , drop=FALSE]
cl   <- cl[match(keep, cl$patient_id), ]

saveRDS(Xt, file.path(PROC,"TCGA_LAML_expression.rds"))
fwrite(cl[,c("sample_id","patient_id","os_time","os_event","age","sex","cohort","platform")],
       file.path(CLIN,"TCGA_LAML_clinical.csv"))
cat("TCGA-LAML: samples with expression =", nrow(Xt), " genes =", ncol(Xt),
    " cases in clinical =", length(cs), "\n")
cat("  vital status table:\n"); print(table(cl$vital, useNA="ifany"))
cat("  usable OS (time>0 & event known):", sum(!is.na(cl$os_time) & cl$os_time>0 & !is.na(cl$os_event)),
    " deaths:", sum(cl$os_event==1 & !is.na(cl$os_time) & cl$os_time>0, na.rm=TRUE), "\n")
