## ============================================================
## 00_setup_dirs.R
## Create the folder skeleton the analysis scripts expect.
## Run once from the repository root:  Rscript scripts/00_setup_dirs.R
## ============================================================
BASE <- Sys.getenv("PROJECT_DIR", unset = getwd())  # repository root
for (d in c("data/raw", "data/processed", "data/clinical",
            "results/tables", "results/models", "results/audit", "logs"))
  dir.create(file.path(BASE, d), recursive = TRUE, showWarnings = FALSE)
cat("Folders ready under", BASE, "\n")
