# Does feature-selection stability predict external reproducibility in AML transcriptomic survival models?

Code, tables and figures for the manuscript *"Does feature-selection stability predict external
reproducibility in AML transcriptomic survival models? A six-cohort benchmark with a permutation
null."* The study uses GSE37642 (GPL570, n = 136) as the discovery cohort and five independent
validation cohorts (GSE37642-GPL96/97, GSE10358, TCGA-LAML, Beat AML, GSE71014).

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

## Data

Raw data are not stored here. Download them into `data/raw/` (created by `00_setup_dirs.R`):

- GSE37642, GSE10358 and GSE71014 from GEO (series matrices, platform annotations, survival file)
- TCGA-LAML STAR counts and file manifest from the GDC
- Beat AML (`aml_ohsu_2022`) from cBioPortal

The exact file names `01_build_cohorts.R` and `02_build_tcga.R` expect are listed at the top of
each script.

## How to run

Run every script from the repository root. Paths are relative to the root, or to the folder set in
the `PROJECT_DIR` environment variable.

```
Rscript scripts/00_setup_dirs.R          # create data/, results/, logs/
Rscript scripts/01_build_cohorts.R       # GEO cohorts -> data/processed, data/clinical
Rscript scripts/02_build_tcga.R          # TCGA-LAML
Rscript scripts/03_audit.R               # cohort audit (Table 1)
Rscript scripts/05_stability.R           # stability vs training size
Rscript scripts/05b_stability_ci.R       # jackknife confidence intervals
Rscript scripts/06_nested_cv.R           # nested cross-validation
Rscript scripts/07_final_external.R      # frozen discovery models, external validation
Rscript scripts/08_sensitivity.R         # sensitivity analyses
Rscript scripts/09_conditional_h4.R      # replication adjusted for discovery effect size
Rscript scripts/10_null_and_swap.R       # permutation null, bootstrap-univariate comparator
Rscript scripts/11_extended_replication.R  # replication in all five validation cohorts
Rscript scripts/12_revalidate_all.R      # external validation, all cohorts (Table 4)
Rscript scripts/13_matched_null.R        # optional: B-matched permutation null
```

`04_hesf_core.R` is the shared function library and is sourced by the other scripts; do not run it
on its own. Scripts 09, 11 and 12 need the outputs of 07, and 05b needs the output of 05. All random
seeds are set in `04_hesf_core.R`. Scripts 05, 06 and 10 are computationally heavy.

## Where each table comes from

Scripts write raw result files to `results/tables/`; the submission tables combine them.

| Submission table | Raw outputs | Script |
|---|---|---|
| Table 1 | `FINAL_COHORT_FEASIBILITY` (cohort audit) | `03_audit.R` (first four cohorts) |
| Table 2 | `TABLE14`, `TABLE20` | `10_null_and_swap.R` |
| Table 3 | `TABLE11`, `TABLE12`, `TABLE15` | `09_conditional_h4.R`, `11_extended_replication.R` |
| Table 4 | `TABLE17` | `12_revalidate_all.R` |
| Table S2 | `TABLE2b`, `TABLE14` | `05_stability.R`, `10_null_and_swap.R` |
| Table S3 | `TABLE9`, `TABLE10` | `08_sensitivity.R` |
| Table S4 | `TABLE3`, nested-CV fold results, `TABLE6` | `06_nested_cv.R`, `07_final_external.R` |
| Table S5 | `TABLE11`, `TABLE12`, `TABLE13`, `TABLE16` | `09_conditional_h4.R`, `11_extended_replication.R` |
| Table S6 | `TABLE4` (top 100 genes) | `07_final_external.R` |
| Table S1 | manual screening of GEO and cBioPortal series | not scripted |

## Notes

- The swapped-discovery analysis in `10_null_and_swap.R` (GSE37642-GPL96/97 as discovery,
  n = 417) was not completed and is not reported.
- The reported permutation null used 25 inner repetitions against 40 in the observed runs. This is
  exact for the single-fit procedures and approximate for the bootstrap-based ones.
- Scripts that build the Beat AML and GSE71014 expression matrices, produce `TABLE16`, `TABLE18` and
  `TABLE20`, draw the figures, and assemble the final submission tables are not yet part of this
  repository.
