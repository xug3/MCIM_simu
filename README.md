# MCIM_simu

This repository contains the R code used for the simulation study accompanying the manuscript:

**“The Missing Covariate Indicator Method is Nearly Valid Almost Always.”**

The simulation study evaluates the performance of the missing covariate indicator method (MCIM) for relative-risk estimation and compares it with:

- analysis of the full data before missingness is imposed (NM);
- multiple imputation (MI); and
- complete-case analysis (CC).

The simulations consider binary exposure, outcome, and covariate variables under missing completely at random (MCAR) and missing at random (MAR) mechanisms. Under MAR, covariate missingness may depend on the exposure, the outcome, or both. The primary simulation outcomes include model convergence, relative bias, and empirical relative efficiency.

## Repository contents

- `batch_LW_10imps.R`  
  Main batch simulation script using 10 imputed datasets for each multiple-imputation analysis.

- `batch_LW_10imps_mar.R`  
  Batch simulation script for the MAR settings.

- `simu_result_countconverg.R`  
  Script for combining simulation results and calculating convergence proportions and summary performance measures.

## Simulation settings

The current implementation uses:

- 2,000 Monte Carlo replicates per parameter setting;
- 10 imputed datasets for each MI analysis;
- log-binomial regression for relative-risk estimation; and
- the R package `jomo` for multiple imputation.

The parameter settings include variations in:

- the prevalence of exposure;
- the prevalence of the partially observed covariate;
- the outcome risk;
- the proportion of covariate missingness;
- the exposure–outcome relative risk;
- the covariate–outcome association;
- the exposure–covariate association; and
- the associations of exposure and outcome with covariate missingness.

## Running the code

1. Download or clone this repository.

2. Open the relevant R script and update any working-directory, input-directory, and output-directory paths.

3. Run the simulation scripts:

   ```r
   source("batch_LW_10imps.R")
   source("batch_LW_10imps_mar.R")
