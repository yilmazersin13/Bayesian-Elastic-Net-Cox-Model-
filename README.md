# Bayesian Elastic-Net Cox Models in Stan (METABRIC Example)

This repository contains a reproducible implementation of a Bayesian elastic-net Cox (BEN–Cox) model for high-dimensional survival data, using the METABRIC breast-cancer cohort as a worked example.

The core idea is to put a global–local elastic-net style prior on the Cox regression coefficients and perform full Bayesian inference with Stan’s HMC sampler. The workflow mirrors the paper:

- preprocess METABRIC data and create a single train/test split,
- fit BEN–Cox on the training set,
- compare against a null Cox model and a ridge-penalised Cox model,
- evaluate IBS, C-index and a Greenwood–Nam–D’Agostino (GND)–style calibration statistic on the held-out test set,
- generate LaTeX-ready tables and PDF figures.

---

## 1. Dependencies

You need a working C++ toolchain for Stan, plus the following R packages:

```r
install.packages(c(
  "survival",
  "ggplot2",
  "dplyr",
  "boot",
  "glmnet",
  "rstan"
))

#EXAMPLE RUN-------------------------------
source("ben_cox_metabric.R")

# quick test run (small iteration counts)
results <- run_complete_analysis(
  data_file = "METABRIC_RNA_Mutation.csv",
  n_iter    = 200,
  n_warmup  = 100
)

# for the paper-level analysis, use something like:
# results <- run_complete_analysis(
#   data_file = "METABRIC_RNA_Mutation.csv",
#   n_iter    = 2000,
#   n_warmup  = 1000
# )
