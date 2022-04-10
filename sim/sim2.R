### Simulation 2: Model Misspecification

simulate2 <- function(scenario, n.sim, a.vals){
  
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
  
  # dr arguments
  bw <- ifelse(n == 800, 0.2, 0.4)
  family <- poisson()
  
  print(scenario)
  
  out <- lapply(1:n.sim, function(i,...){
    
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
    
    # naive
    naive_hat <- try(erf(y = y, a = z.tilde, x = x, offset = offset, weights = weights, 
                         family = family, a.vals = a.vals, bw = bw,
                         n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # real
    rc_hat <- try(erf(y = y, a = z.hat, x = x, offset = offset, weights = weights,
                      family = family, a.vals = a.vals, bw = bw,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # BART Approach
    bart_hat <- try(bart_erf(s = s, s.tilde = s.tilde, y = y, s.id = s.id, id = id, 
                             w = w, x = x, offset = offset, a.vals = a.vals, bw = bw,
                             scale = scale, shape = shape, rate = rate, 
                             n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # Bayes DR Approach
    bayes_hat <- try(bayes_erf(s = s, s.tilde = s.tilde, y = y, s.id = s.id, id = id,
                               w = w, x = x, offset = offset, a.vals = a.vals, bw = bw,
                               scale = scale, shape = shape, rate = rate, dr = TRUE,
                               n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # estimates
    est <- rbind(predict_example(a = a.vals, x = x, out_scen = out_scen),
                 if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(rc_hat, "try-error")) {rc_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bart_hat, "try-error")) {bart_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bayes_hat, "try-error")) {bayes_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bayes_hat, "try-error")) {sqrt(bayes_hat$estimate.dr)} else {rep(NA, length(a.vals))})
    
    #standard error
    se <- rbind(if (!inherits(naive_hat, "try-error")) {sqrt(naive_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(rc_hat, "try-error")) {sqrt(rc_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bayes_hat, "try-error")) {sqrt(bayes_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bayes_hat, "try-error")) {sqrt(bayes_hat$variance.dr)} else {rep(NA, length(a.vals))})
    
    return(list(est = est, se = se))
    
  }s)
  
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
  rownames(out_est) <- c("ERF","NAIVE","RC","BART","GLM","DR")
  
  out_bias <- t(apply(est[2:6,,], 1, function(x) rowMeans(abs(x - mu.mat), na.rm = T)))
  colnames(out_bias) <- a.vals
  rownames(out_bias) <- c("NAIVE","RC","BART","GLM","DR")
  
  out_mse <- t(apply(est[2:6,,], 1, function(x) rowMeans((x - mu.mat)^2, na.rm = T)))
  colnames(out_mse) <- a.vals
  rownames(out_mse) <- c("NAIVE","RC","BART","GLM","DR")
  
  out_cp <- do.call(rbind, lapply(cp, rowMeans, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("NAIVE","RC","BART","GLM","DR")
  
  lower <- rbind(rowMeans(est[2,,] - 1.96*se[1,,]),
                 rowMeans(est[3,,] - 1.96*se[2,,]),
                 rowMeans(est[4,,] - 1.96*se[3,,]),
                 rowMeans(est[5,,] - 1.96*se[4,,]),
                 rowMeans(est[6,,] - 1.96*se[5,,]))
  
  colnames(lower) <- a.vals
  rownames(lower) <- c("NAIVE","RC","BART","GLM","DR")
  
  upper <- rbind(rowMeans(est[2,,] + 1.96*se[1,,]),
                 rowMeans(est[3,,] + 1.96*se[2,,]),
                 rowMeans(est[4,,] + 1.96*se[3,,]),
                 rowMeans(est[5,,] + 1.96*se[4,,]),
                 rowMeans(est[6,,] + 1.96*se[5,,]))
                 
  colnames(upper) <- a.vals
  rownames(upper) <- c("NAIVE","RC","BART","GLM","DR")

  
  rslt <- list(scenario = scenario, est = out_est, bias = out_bias, mse = out_mse, cp = out_cp, lower = lower, upper = upper)
  filename <- paste0("~/Dropbox/Projects/ERC-EPE/Output/sim_2/", paste(scenario, collapse = "_"),".RData")
  save(rslt, file = filename)
  
}

rm(list = ls())

## Preliminaries

library(parallel)

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 14, by = 0.08)
n.sim <- 500

n <- 800
mult <- 5
sig_agg <- sqrt(2)
sig_pred <- 1
gps_scen <- c("a", "b")
out_scen <- c("a", "b")
pred_scen <- c("a", "b")

scen_mat <- expand.grid(n = n, mult = mult, sig_agg = sig_agg, sig_pred = sig_pred, 
                        gps_scen = gps_scen, out_scen = out_scen, pred_scen = pred_scen, 
                        stringsAsFactors = FALSE)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- mclapply(scenarios, simulate2, n.sim = n.sim, a.vals = a.vals, mc.cores = 8, mc.preschedule = TRUE)
