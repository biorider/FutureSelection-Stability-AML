## ============================================================
## HESF-AML  |  01_build_cohorts.R
## Parse GEO series matrices + TCGA-LAML into analysis-ready
## expression matrices (patients x genes) and clinical tables.
## ============================================================
suppressPackageStartupMessages({library(data.table)})

BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
RAW  <- file.path(BASE,"data","raw");       PROC <- file.path(BASE,"data","processed")
CLIN <- file.path(BASE,"data","clinical");  AUD  <- file.path(BASE,"results","audit")
for (d in c(PROC,CLIN,AUD)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

## ---------- helpers -------------------------------------------------
gunzip_to <- function(gz, out){
  if (file.exists(out) && file.size(out) > file.size(gz)) return(out)
  ic <- gzfile(gz,"rt"); oc <- file(out,"wt")
  repeat { l <- readLines(ic, n=50000, warn=FALSE); if(!length(l)) break; writeLines(l, oc) }
  close(ic); close(oc); out
}

read_series <- function(gz){
  txt <- gunzip_to(gz, file.path(tempdir(), sub("\\.gz$","",basename(gz))))
  L <- readLines(txt, n=500, warn=FALSE)
  b <- grep("^!series_matrix_table_begin", L)[1]
  hdr <- L[1:(b-1)]
  X <- fread(txt, skip=b, header=TRUE, fill=TRUE, showProgress=FALSE,
             na.strings=c("NA","null",""), data.table=FALSE)
  X <- X[!startsWith(as.character(X[[1]]), "!"), , drop=FALSE]
  ids <- as.character(X[[1]]); X <- X[,-1,drop=FALSE]
  M <- as.matrix(X); storage.mode(M) <- "double"; rownames(M) <- ids
  colnames(M) <- gsub('"','',colnames(M))
  list(hdr=hdr, X=M)
}

meta_of <- function(hdr){
  g1 <- function(tag) gsub('"','',strsplit(grep(tag,hdr,value=TRUE)[1],"\t")[[1]])[-1]
  rows <- grep("^!Sample_characteristics_ch1", hdr, value=TRUE)
  ch <- do.call(rbind, lapply(rows, function(r) gsub('"','',strsplit(r,"\t")[[1]])[-1]))
  list(gsm=g1("^!Sample_geo_accession"), title=g1("^!Sample_title"),
       platform=g1("^!Sample_platform_id"), ch=ch)
}
getch <- function(ch, key){
  k <- paste0(tolower(key),":")
  apply(ch, 2, function(col){ v <- col[startsWith(tolower(col), k)]
    if (length(v)) trimws(sub("^[^:]*:","",v[1])) else NA_character_ })
}

## GEO platform annotation -> probe/gene symbol map
read_gpl <- function(file){
  L <- readLines(file, warn=FALSE)
  b <- grep("^!platform_table_begin", L)[1]
  e <- grep("^!platform_table_end",   L)[1]
  tb <- fread(text=L[(b+1):(e-1)], sep="\t", header=TRUE, quote="",
              showProgress=FALSE, data.table=FALSE)
  sym <- grep("^Gene Symbol$|^Gene symbol$|^GENE_SYMBOL$", names(tb), value=TRUE)[1]
  if (is.na(sym)) sym <- grep("symbol", names(tb), ignore.case=TRUE, value=TRUE)[1]
  out <- data.frame(probe=as.character(tb[[1]]), symbol=trimws(as.character(tb[[sym]])),
                    stringsAsFactors=FALSE)
  ## probesets mapping to several genes ("A /// B") are ambiguous -> dropped
  out$symbol[grepl("///", out$symbol, fixed=TRUE)] <- NA
  out$symbol[out$symbol %in% c("","---","NA")] <- NA
  out[!is.na(out$symbol), ]
}

## probe matrix (probes x samples) -> gene matrix (samples x genes)
## prespecified aggregation: per gene keep probeset with highest mean
## expression across that cohort (unsupervised; no outcome involved)
collapse_to_genes <- function(M, map){
  map <- map[map$probe %in% rownames(M), ]
  mu  <- rowMeans(M[map$probe, , drop=FALSE], na.rm=TRUE)
  ord <- order(map$symbol, -mu)
  map <- map[ord, ]
  keep <- map[!duplicated(map$symbol), ]
  G <- t(M[keep$probe, , drop=FALSE])
  colnames(G) <- keep$symbol
  G
}

log2_if_needed <- function(M, label){
  rng <- range(M, na.rm=TRUE)
  if (rng[2] > 100) { M[M < 1] <- 1; M <- log2(M)
    message(label, ": linear scale detected (max=", round(rng[2],1), ") -> log2 transformed") }
  else message(label, ": already log-like (range ", paste(round(rng,2),collapse=" .. "), ")")
  M
}

## ---------- GPL annotation -------------------------------------------
map570 <- read_gpl(file.path(RAW,"GPL570_annot.txt"))
map96  <- read_gpl(file.path(RAW,"GPL96_annot.txt"))
map97  <- read_gpl(file.path(RAW,"GPL97_annot.txt"))
saveRDS(list(GPL570=map570, GPL96=map96, GPL97=map97), file.path(PROC,"probe_maps.rds"))

## ---------- Cohort 1: GSE37642 / GPL570 (discovery) -------------------
s570 <- read_series(file.path(RAW,"GSE37642-GPL570_series_matrix.txt.gz"))
m570 <- meta_of(s570$hdr)
cl570 <- data.frame(
  sample_id  = m570$gsm,
  patient_id = m570$title,
  os_time    = suppressWarnings(as.numeric(getch(m570$ch,"overall survival (days)"))),
  os_event   = ifelse(tolower(getch(m570$ch,"life status"))=="dead",1L,
               ifelse(tolower(getch(m570$ch,"life status"))=="alive",0L,NA_integer_)),
  age        = suppressWarnings(as.numeric(getch(m570$ch,"age"))),
  fab        = getch(m570$ch,"fab"),
  runx1_runx1t1_fusion = getch(m570$ch,"runx1-runx1t1_fusion"),
  runx1_mutation       = getch(m570$ch,"runx1_mutation"),
  platform   = m570$platform, cohort = "GSE37642_GPL570",
  stringsAsFactors = FALSE)

## cross-check against the GEO supplementary survival file
sv <- fread(gunzip_to(file.path(RAW,"GSE37642_Survival_data.txt.gz"),
            file.path(tempdir(),"GSE37642_Survival_data.txt")),
            sep="\t", header=TRUE, data.table=FALSE)
names(sv) <- c("sample_id","patient_id","os_days_suppl","status_suppl")
chk <- merge(cl570, sv, by="sample_id", all.x=TRUE)
sup_check <- data.frame(
  n_matrix_samples   = nrow(cl570),
  n_in_suppl_file    = sum(!is.na(chk$os_days_suppl)),
  n_os_time_agree    = sum(chk$os_time == chk$os_days_suppl, na.rm=TRUE),
  n_patient_id_agree = sum(chk$patient_id.x == chk$patient_id.y, na.rm=TRUE),
  n_suppl_total_rows = nrow(sv))
fwrite(sup_check, file.path(AUD,"GSE37642_supplementary_crosscheck.csv"))

## ---------- Cohort 2: GSE37642 / GPL96+GPL97 (platform-shift) ---------
s96 <- read_series(file.path(RAW,"GSE37642-GPL96_series_matrix.txt.gz"))
s97 <- read_series(file.path(RAW,"GSE37642-GPL97_series_matrix.txt.gz"))
m96 <- meta_of(s96$hdr); m97 <- meta_of(s97$hdr)
cl96 <- data.frame(
  sample_id  = m96$gsm, patient_id = m96$title,
  os_time    = suppressWarnings(as.numeric(getch(m96$ch,"overall survival (days)"))),
  os_event   = ifelse(tolower(getch(m96$ch,"life status"))=="dead",1L,
               ifelse(tolower(getch(m96$ch,"life status"))=="alive",0L,NA_integer_)),
  age        = suppressWarnings(as.numeric(getch(m96$ch,"age"))),
  fab        = getch(m96$ch,"fab"),
  platform   = "GPL96_GPL97", cohort = "GSE37642_GPL96_97",
  stringsAsFactors = FALSE)

## ---------- Cohort 3: GSE10358 (population-shift) ---------------------
s103 <- read_series(file.path(RAW,"GSE10358_series_matrix.txt.gz"))
m103 <- meta_of(s103$hdr)
os_m <- suppressWarnings(as.numeric(getch(m103$ch,"os months (3.31.10)")))
vs   <- getch(m103$ch,"vital status")
cl103 <- data.frame(
  sample_id  = m103$gsm, patient_id = m103$gsm,
  os_time    = os_m * 30.4375,
  os_event   = ifelse(grepl("^dead", tolower(vs)),1L,
               ifelse(grepl("^alive",tolower(vs)),0L,NA_integer_)),
  age        = suppressWarnings(as.numeric(getch(m103$ch,"age on study"))),
  sex        = getch(m103$ch,"Sex"),
  fab        = getch(m103$ch,"fab"),
  platform   = m103$platform, cohort = "GSE10358",
  stringsAsFactors = FALSE)

## ---------- raw probe-level audit before collapsing -------------------
plat_summary <- data.frame(
  dataset  = c("GSE37642-GPL570","GSE37642-GPL96","GSE37642-GPL97","GSE10358"),
  platform = c("GPL570","GPL96","GPL97","GPL570"),
  n_samples  = c(ncol(s570$X), ncol(s96$X), ncol(s97$X), ncol(s103$X)),
  n_features = c(nrow(s570$X), nrow(s96$X), nrow(s97$X), nrow(s103$X)),
  missing_values = c(sum(is.na(s570$X)), sum(is.na(s96$X)), sum(is.na(s97$X)), sum(is.na(s103$X))),
  min_value = round(c(min(s570$X,na.rm=TRUE),min(s96$X,na.rm=TRUE),min(s97$X,na.rm=TRUE),min(s103$X,na.rm=TRUE)),3),
  max_value = round(c(max(s570$X,na.rm=TRUE),max(s96$X,na.rm=TRUE),max(s97$X,na.rm=TRUE),max(s103$X,na.rm=TRUE)),3))
fwrite(plat_summary, file.path(AUD,"GEO_platform_summary.csv"))

## ---------- collapse to genes ----------------------------------------
X570 <- collapse_to_genes(log2_if_needed(s570$X,"GSE37642-GPL570"), map570)
X103 <- collapse_to_genes(log2_if_needed(s103$X,"GSE10358"),        map570)
## GPL96 + GPL97 = HG-U133A + HG-U133B for the SAME patients
stopifnot(identical(m96$title, m97$title))
X96 <- collapse_to_genes(log2_if_needed(s96$X,"GSE37642-GPL96"), map96)
X97 <- collapse_to_genes(log2_if_needed(s97$X,"GSE37642-GPL97"), map97)
rownames(X96) <- m96$title; rownames(X97) <- m97$title
common_ab <- setdiff(colnames(X97), colnames(X96))     # B-chip adds new genes only
X9697 <- cbind(X96, X97[rownames(X96), common_ab, drop=FALSE])
rownames(X570) <- cl570$patient_id
rownames(X103) <- cl103$patient_id
rownames(cl96) <- NULL

saveRDS(X570,  file.path(PROC,"GSE37642_GPL570_expression.rds"))
saveRDS(X9697, file.path(PROC,"GSE37642_GPL96_97_expression.rds"))
saveRDS(X103,  file.path(PROC,"GSE10358_expression.rds"))
fwrite(cl570, file.path(CLIN,"GSE37642_GPL570_clinical.csv"))
fwrite(cl96,  file.path(CLIN,"GSE37642_GPL96_97_clinical.csv"))
fwrite(cl103, file.path(CLIN,"GSE10358_clinical.csv"))

cat("\n--- built ---\n")
print(plat_summary)
print(sup_check)
cat("genes: GPL570=",ncol(X570)," GPL96_97=",ncol(X9697)," GSE10358=",ncol(X103),"\n",sep="")
