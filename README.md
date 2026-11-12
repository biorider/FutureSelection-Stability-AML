# Does feature-selection stability predict external reproducibility in AML transcriptomic survival models?

Code, tables and figures for the manuscript *"Does feature-selection stability predict external
reproducibility in AML transcriptomic survival models?"* The study uses GSE37642 (GPL570, n = 136) as the discovery cohort and five independent validation cohorts (GSE37642-GPL96/97, GSE10358, TCGA-LAML, Beat AML, GSE71014).

## Repository contents

| Folder | Contents |
|---|---|
| `scripts/` | R scripts that build the cohorts and run every analysis (see run order below) |
| `tables/` | Main Tables 1–4 (`Tables_1-4.docx`) and Supplementary Tables S1–S6 (one `.xlsx` each) |
| `figures/` | Main Figures 1–7 and Supplementary Figure S1 (PNG) |

## Requirements

R with the packages `survival`, `glmnet`, `ranger`, `matrixStats`, `data.table`, `jsonlite`,
`timeROC` and `pROC`. The stability experiments run on a parallel cluster; several scripts call
`make_hesf_cluster(18)`, so change `18` to the number of cores on your machine.

