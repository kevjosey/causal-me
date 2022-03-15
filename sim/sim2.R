### Test DR estimator

rm(list = ls())

## Preliminaries
library(mvtnorm)
library(SuperLearner)
library(parallel)
library(abind)
library(dbarts)

# Code for generating and fitting data
source("~/Github/causal-me/sim/gen-data.R")
source("~/Github/causal-me/erf.R")
source("~/Github/causal-me/erf-alt.R")
source("~/Github/causal-me/bart-erf.R")
source("~/Github/causal-me/bart-erf-alt.R")
source("~/Github/causal-me/bayes-erf.R")
source("~/Github/causal-me/auxiliary.R")

simulate <- function(scenario, n.sim, a.vals){
  
  # simulation arguments
  sig_gps <- 2
  sig_agg <- scenario$sig_agg
  sig_pred <- scenario$sig_pred
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  pred_scen <- scenario$pred_scen
  prob <- 0.1
  
  # gen data arguments
  n <- scenario$n # c(500, 800)
  mult <- scenario$mult # c(100, 200)
  
  # gibbs sampler stuff
  n.iter <- 2000
  n.adapt <- 2000
  thin <- 20
  h.a <- 0.5
  scale <- 1e6
  shape <- rate <- 1e-3
  
  # dr arguments
  span <- ifelse(n == 800, 0.125, 0.25)
  family <- poisson()
  
  print(scenario)
  
  out <- mclapply(1:n.sim, function(i,...){
    
    dat <- gen_data(n = n, mult = mult, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
                    pred_scen = pred_scen, out_scen = out_scen, gps_scen = gps_scen)
    
    # zipcode index
    s.id <- dat$s.id
    id <- dat$id
    offset <- log(dat$offset)
    weights <- dat$weights
    
    # data
    y <- dat$y
    x <- dat$x
    a <- dat$a
    w <- dat$w
    s_tilde <- dat$star
    
    # validation subset
    s <- dat$s*rbinom(mult*n, 1, prob)
    s[s == 0] <- NA
    
    # remove clusters w/o exposure data
    keep <- which((id %in% unique(s.id)))
    if (length(keep)!= 0) {
      y <- dat$y[(id %in% keep)]
      x <- dat$x[(id %in% keep),]
      a <- dat$a[(id %in% keep)]
      offset <- log(dat$offset[(id %in% keep)])
      id <- dat$id[(id %in% keep)]
    }
    
    # exposure predictions
    s_hat <- pred(s = s, star = s_tilde, w = w, sl.lib = "SL.glm")
    z_hat <- aggregate(s_hat, by = list(s.id), mean)[,2]
    z_tilde <- aggregate(s_tilde, by = list(s.id), mean)[,2]
    
    # naive
    naive_hat <- try(erf(y = y, a = z_tilde, x = x, offset = offset, weights = weights, 
                         family = family, a.vals = a.vals, span = span,
                         n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # real
    rc_hat <- try(erf(y = y, a = z_hat, x = x, offset = offset, weights = weights,
                      family = family, a.vals = a.vals, span = span,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # BART Approach
    bart_hat <- try(bart_erf(s = s, star = s_tilde, y = y, offset = offset, weights = weights,
                             s.id = s.id, id = id, w = w, x = x, family = family,
                             a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                             h.a = h.a, n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # Alternate Bayes
    bayes_hat <- try(bayes_erf(s = s, star = s_tilde, y = y, offset = offset, weights = weights,
                               s.id = s.id, id = id, w = w, x = x, family = family, df = 6,
                               a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                               h.a = h.a, n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # estimates
    est <- rbind(predict_example(a = a.vals, x = x, out_scen = out_scen),
                 if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(rc_hat, "try-error")) {rc_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bart_hat, "try-error")) {bart_hat$smooth_estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bayes_hat, "try-error")) {bayes_hat$dr_estimate} else {rep(NA, length(a.vals))})
    
    #standard error
    se <- rbind(if (!inherits(naive_hat, "try-error")) {sqrt(naive_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(rc_hat, "try-error")) {sqrt(rc_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$smooth_variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bayes_hat, "try-error")) {sqrt(bayes_hat$dr_variance)} else {rep(NA, length(a.vals))})
    
    return(list(est = est, se = se))
    
  }, mc.cores = 30, mc.preschedule = TRUE)
  
  est <- abind(lapply(out, function(lst, ...) if (!inherits(lst, "try-error")) {lst$est} else {matrix(NA, ncol = length(a.vals), nrow = 6)}), along = 3)
  se <- abind(lapply(out, function(lst, ...) if (!inherits(lst, "try-error")) {lst$se} else {matrix(NA, ncol = length(a.vals), nrow = 5)}), along = 3)
  mu.mat <- matrix(rep(rowMeans(est[1,,]), n.sim), nrow = length(a.vals), ncol = n.sim)
  
  # coverage probability
  cp <- list(as.matrix((est[2,,] - 1.96*se[1,,]) < mu.mat & (est[2,,] + 1.96*se[1,,]) > mu.mat),
             as.matrix((est[3,,] - 1.96*se[2,,]) < mu.mat & (est[3,,] + 1.96*se[2,,]) > mu.mat),
             as.matrix((est[4,,] - 1.96*se[3,,]) < mu.mat & (est[4,,] + 1.96*se[3,,]) > mu.mat),
             as.matrix((est[5,,] - 1.96*se[4,,]) < mu.mat & (est[5,,] + 1.96*se[4,,]) > mu.mat),
             as.matrix((est[6,,] - 1.96*se[5,,]) < mu.mat & (est[6,,] + 1.96*se[5,,]) > mu.mat))
  
  out_est <- t(apply(est, 1, rowMeans, na.rm = T))
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERF","NAIVE","RC","BART","LOESS","DR")
  
  out_bias <- t(apply(est[2:6,,], 1, function(x) rowMeans(abs(x - mu.mat), na.rm = T)))
  colnames(out_bias) <- a.vals
  rownames(out_bias) <- c("NAIVE","RC","BART","LOESS","DR")
  
  out_mse <- t(apply(est[2:6,,], 1, function(x) rowMeans((x - mu.mat)^2, na.rm = T)))
  colnames(out_mse) <- a.vals
  rownames(out_mse) <- c("NAIVE","RC","BART","LOESS","DR")
  
  out_cp <- do.call(rbind, lapply(cp, rowMeans, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("NAIVE","RC","BART","LOESS","DR")
  
  lower <- rbind(rowMeans(est[2,,] - 1.96*se[1,,]),
                 rowMeans(est[3,,] - 1.96*se[2,,]),
                 rowMeans(est[4,,] - 1.96*se[3,,]),
                 rowMeans(est[5,,] - 1.96*se[4,,]),
                 rowMeans(est[6,,] - 1.96*se[5,,]))
  
  colnames(lower) <- a.vals
  rownames(lower) <- c("NAIVE","RC","BART","LOESS","DR")
  
  upper <- rbind(rowMeans(est[2,,] + 1.96*se[1,,]),
                 rowMeans(est[3,,] + 1.96*se[2,,]),
                 rowMeans(est[4,,] + 1.96*se[3,,]),
                 rowMeans(est[5,,] + 1.96*se[4,,]),
                 rowMeans(est[6,,] + 1.96*se[5,,]))
                 
  colnames(upper) <- a.vals
  rownames(upper) <- c("NAIVE","RC","BART","LOESS","DR")

  
  rslt <- list(scenario = scenario, est = out_est, bias = out_bias, mse = out_mse, cp = out_cp, lower = lower, upper = upper)
  filename <- paste0("~/Dropbox/Projects/ERF-EPE/Output/sim_2/", paste(scenario, collapse = "_"),".RData")
  save(rslt, file = filename)
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 14, by = 0.04)
n.sim <- 500

n <- 800
mult <- 5
sig_agg <- sqrt(2)
sig_pred <- 1
gps_scen <- c("a", "b")
out_scen <- c("a", "b")
pred_scen <- "b"

scen_mat <- expand.grid(n = n, mult = mult, sig_agg = sig_agg, sig_pred = sig_pred, gps_scen = gps_scen, out_scen = out_scen, pred_scen = pred_scen, stringsAsFactors = FALSE)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- lapply(scenarios, simulate, n.sim = n.sim, a.vals = a.vals)
