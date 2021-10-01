# Estimating a Causal Exposure Response Function with a Continuous Error-Prone Exposure: A Study of Fine Particulate Matter and All-Cause Mortality 

[`erf.R`](https://github.com/kevjosey/causal-me/tree/master/erf.R) Includes baseline functions for fitting an exposure response function (ERF) without measurement error correction. Code in this script is later used by `bart-erf.R`.

[`erf-alt.R`](https://github.com/kevjosey/causal-me/tree/master/erf-alt.R) An alternative implementation to `erf.R` where in addition the exposure generalized propensity score is predicted and the ERF estimator incorporates doubly-robust pseudo-outcomes. These functions are later used by `bayes-erf.R` and `bayes-erf-alt`.

[`bart-erf.R`](https://github.com/kevjosey/causal-me/tree/master/bart-erf.R) Includes the function `bart_erf()` which fits a measurement error corrected ERF using multiple imputation and a Bayesian additive regression tree (BART) outcome model.

[`bart-erf-alt.R`](https://github.com/kevjosey/causal-me/tree/master/bart-erf-alt.R) An alternative function to `bart-erf.R` which instead projects the doubly-robust pseudo-outcome onto the support of the exposure.

[`bayes-erf.R`](https://github.com/kevjosey/causal-me/tree/master/bayes-erf.R) An alternative function to `bart-erf-alt.R` which fits a generalized linear outcome model before regressing the pseudo-outcome onto the exposure support.

[`auxiliary.R`](https://github.com/kevjosey/causal-me/tree/master/auxiliary.R) Additional functions used intermittently throughout the manuscript including a function to estimate the regression calibrated exposures and a function to compute the highest posterior density interval when the BART ERF is used without smoothing.

[`/sim`](https://github.com/kevjosey/causal-me/tree/mastergen-data.R) Code for running the numerical studies is contained in this directory. Also included in this folder is `test.R` for unit testing and `gen-data.R` which includes functions to generate simulated data and the true exposure response function.
