# Estimating a Causal Exposure Response Function with a Continuous Error-Prone Exposure: A Study of Fine Particulate Matter and All-Cause Mortality 

[`erf.R`](https://github.com/kevjosey/causal-me/tree/master/erf.R) Includes baseline functions for fitting an exposure response function (ERF) without measurement error correction. Code in this script is later used by `bart-erf.R`, `bayes-erf.R`, and `spatial-erf.R`.

[`bart-erf.R`](https://github.com/kevjosey/causal-me/tree/master/bart-erf.R) Includes the function `bart_erf()` which fits a measurement error corrected ERF using multiple imputation and a Bayesian additive regression tree (BART) outcome model.

[`bayes-erf.R`](https://github.com/kevjosey/causal-me/tree/master/bayes-erf.R) An alternative function to `bart_erf()` which fits a generalized linear  model of the outcome before regressing the pseudo-outcome onto the exposure support.

[`spatial-erf.R`](https://github.com/kevjosey/causal-me/tree/master/spatial-erf.R) An alternative function to `bart_erf()` which incorporates a spatial autocorrelation random variable into the cluster-level exposure model.

[`auxiliary.R`](https://github.com/kevjosey/causal-me/tree/master/auxiliary.R) Additional functions used intermittently throughout the manuscript including a function to estimate the regression calibrated grid-level exposures, a function to compute the highest posterior density interval when the `bart_erf()` and `bayes_erf()` is used without smoothing, and a function to check the adjacency matrix for the spatial autocorrelation component.

[`/sim`](https://github.com/kevjosey/causal-me/tree/master/sim) Code for running the numerical studies is contained in this directory. Also included in this folder is `test.R` for unit testing and `gen-data.R` which includes functions to generate simulated data and the true exposure response function.
