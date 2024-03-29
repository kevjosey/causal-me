### Simulation 1: Measurement Error

## Preliminaries

rm(list = ls())

library(parallel)
library(mvtnorm)
library(abind)
library(SuperLearner)
library(dbarts)

# Code for generating and fitting data
source("~/Github/causal-me/sim/gen-data.R")
source("~/Github/causal-me/erf.R")
source("~/Github/causal-me/bart-erf.R")
source("~/Github/causal-me/bayes-erf.R")
source("~/Github/causal-me/auxiliary.R")

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 14, by = 0.04)
n.sim <- 500

n <- c(400, 800)
mult <- c(5, 10)
sig_agg <- c(0, 1, sqrt(2))
sig_pred <- c(0, 1, sqrt(2))
gps_scen <- "a"
out_scen <- "a"
pred_scen <- "a"

scen_mat <- expand.grid(n = n, mult = mult, sig_agg = sig_agg, sig_pred = sig_pred,
                        gps_scen = gps_scen, out_scen = out_scen, pred_scen = pred_scen, 
                        stringsAsFactors = FALSE)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])

for (i in 1:length(scenarios)) {
  
  scenario <- scenarios[[i]]
  
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
  scale <- 1e6
  shape <- 1e-3
  rate <- 1e-3
  bw <- ifelse(n == 800, 0.2, 0.4)
  family <- poisson()
  
  print(scenario)
  
  out <- mclapply(1:n.sim, function(i,...){
    
    dat <- gen_data(n = n, mult = mult, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
                    pred_scen = pred_scen, out_scen = out_scen, gps_scen = gps_scen)
    
    # zipcode index
    s.id <- dat$s.id
    id <- dat$id
    offset <- dat$offset
    weights <- dat$weights
    
    # data
    y <- dat$y
    x <- dat$x
    a <- dat$a
    w <- dat$w
    s.tilde <- dat$s.tilde
    
    # validation subset
    s <- dat$s*rbinom(mult*n, 1, prob)
    s[s == 0] <- NA
    
    # remove clusters w/o exposure data
    keep <- which((id %in% unique(s.id)))
    if (length(keep)!= 0) {
      y <- dat$y[(id %in% keep)]
      x <- dat$x[(id %in% keep),]
      a <- dat$a[(id %in% keep)]
      offset <- dat$offset[(id %in% keep)]
      id <- dat$id[(id %in% keep)]
    }
    
    # exposure predictions
    s.hat <- pred(s = s, s.tilde = s.tilde, w = w, sl.lib = "SL.glm")
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    z.tilde <- aggregate(s.tilde, by = list(s.id), mean)[,2]
    
    # real
    rc_hat <- try(erf(y = y, a = z.hat, x = x, offset = offset, a.vals = a.vals, bw = bw), silent = TRUE)
    
    # naive
    naive_hat <- try(erf(y = y, a = z.tilde, x = x, offset = offset, a.vals = a.vals, bw = bw), silent = TRUE)
    
    # BART Approach
    bart_hat <- try(bart_erf(s = s, s.tilde = s.tilde, y = y, s.id = s.id, id = id, 
                             w = w, x = x, offset = offset, a.vals = a.vals, bw = bw,
                             scale = scale, shape = shape, rate = rate, 
                             n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # Bayes DR Approach
    bayes_hat <- try(bayes_erf(s = s, s.tilde = s.tilde, y = y, s.id = s.id, id = id,
                               w = w, x = x, offset = offset, a.vals = a.vals, bw = bw,
                               scale = scale, shape = shape, rate = rate,
                               n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # estimates
    est <- rbind(predict_example(a = a.vals, x = x, offset = offset, out_scen = out_scen),
                 if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(rc_hat, "try-error")) {rc_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bart_hat, "try-error")) {bart_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bayes_hat, "try-error")) {bayes_hat$estimate} else {rep(NA, length(a.vals))})
    
    #standard error
    se <- rbind(if (!inherits(naive_hat, "try-error")) {sqrt(naive_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(rc_hat, "try-error")) {sqrt(rc_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bayes_hat, "try-error")) {sqrt(bayes_hat$variance)} else {rep(NA, length(a.vals))})
    
    return(list(est = est, se = se))
    
  }, mc.cores = 30)
  
  est <- abind(lapply(out, function(lst, ...) if (!inherits(lst, "try-error")) {lst$est} else {matrix(NA, ncol = length(a.vals), nrow = 5)}), along = 3)
  se <- abind(lapply(out, function(lst, ...) if (!inherits(lst, "try-error")) {lst$se} else {matrix(NA, ncol = length(a.vals), nrow = 4)}), along = 3)
  mu.mat <- matrix(rep(rowMeans(est[1,,]), n.sim), nrow = length(a.vals), ncol = n.sim)
  
  # coverage probability
  cp <- list(as.matrix((est[2,,] - 1.96*se[1,,]) < mu.mat & (est[2,,] + 1.96*se[1,,]) > mu.mat),
             as.matrix((est[3,,] - 1.96*se[2,,]) < mu.mat & (est[3,,] + 1.96*se[2,,]) > mu.mat),
             as.matrix((est[4,,] - 1.96*se[3,,]) < mu.mat & (est[4,,] + 1.96*se[3,,]) > mu.mat),
             as.matrix((est[5,,] - 1.96*se[4,,]) < mu.mat & (est[5,,] + 1.96*se[4,,]) > mu.mat))
  
  out_est <- t(apply(est, 1, rowMeans, na.rm = T))
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERF","NAIVE","RC","BART","GLM")
  
  out_bias <- t(apply(est[2:5,,], 1, function(x) rowMeans(abs(x - mu.mat), na.rm = T)))
  colnames(out_bias) <- a.vals
  rownames(out_bias) <- c("NAIVE","RC","BART","GLM")
  
  out_mse <- t(apply(est[2:5,,], 1, function(x) rowMeans((x - mu.mat)^2, na.rm = T)))
  colnames(out_mse) <- a.vals
  rownames(out_mse) <- c("NAIVE","RC","BART","GLM")
  
  out_cp <- do.call(rbind, lapply(cp, rowMeans, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("NAIVE","RC","BART","GLM")
  
  lower <- rbind(rowMeans(est[2,,] - 1.96*se[1,,]),
                 rowMeans(est[3,,] - 1.96*se[2,,]),
                 rowMeans(est[4,,] - 1.96*se[3,,]),
                 rowMeans(est[5,,] - 1.96*se[4,,]))
  
  colnames(lower) <- a.vals
  rownames(lower) <- c("NAIVE","RC","BART","GLM")
  
  upper <- rbind(rowMeans(est[2,,] + 1.96*se[1,,]),
                 rowMeans(est[3,,] + 1.96*se[2,,]),
                 rowMeans(est[4,,] + 1.96*se[3,,]),
                 rowMeans(est[5,,] + 1.96*se[4,,]))
  
  colnames(upper) <- a.vals
  rownames(upper) <- c("NAIVE","RC","BART","GLM")
  
  rslt <- list(scenario = scenario, est = out_est, bias = out_bias, mse = out_mse, cp = out_cp, lower = lower, upper = upper)
  filename <- paste0("~/Dropbox/Projects/ERC-EPE/Output/sim_1/", paste(scenario, collapse = "_"),".RData")
  save(rslt, file = filename)
  
}
