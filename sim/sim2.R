### Test DR estimator

rm(list = ls())

## Preliminaries

library(mvtnorm)
library(SuperLearner)
library(parallel)
library(abind)
library(dbarts)

# Code for generating and fitting data
source("~/Github/causal-me/auxiliary.R")
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/bart-erc.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/erc.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # simulation arguments
  sig_gps <- scenario$sig_gps
  sig_agg <- sqrt(2)
  sig_pred <- sqrt(0.5)
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  pred_scen <- scenario$pred_scen
  prob <- 0.1
  
  # gen data arguments
  n <- scenario$n # c(500, 800)
  mult <- scenario$mult # c(100, 200)
  
  # gibbs sampler stuff
  thin <- 20
  n.iter <- 10000
  n.adapt <- 1000
  h.a <- 1
  h.gamma <- 0.03
  scale <- 1e6
  shape <- rate <- 1e-3
  deg.num <- 2
  
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
    s_hat <- pred(s = s, star = s_tilde, w = w, sl.lib = sl.lib)
    a_hat <- blp(s = s_hat, s.id = s.id, x = x, id = id)$a
    z_hat <- aggregate(s_hat, by = list(s.id), mean)[,2]
    
    # real
    obs_hat <- try(erc(y = y, a = a, x = x, offset = offset, family = family,
                       a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    
    # naive
    naive_hat <- try(erc(y = y, a = z_hat, x = x, offset = offset, family = family,
                         a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num))
    
    # Bayesian Approach
    bayes_hat <- try(glm_erc(s = s, star = s_tilde, y = y, offset = offset,
                             s.id = s.id, id = id, w = w, x = x, family = family,
                             a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                             h.a = h.a, h.gamma = h.gamma, n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # BART Approach
    bart_hat <- try(bart_erc(s = s, star = s_tilde, y = y, offset = offset,
                             s.id = s.id, id = id, w = w, x = x, family = family,
                             a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                             h.a = h.a, n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # estimates
    est <- rbind(predict_example(a = a.vals, x = x, out_scen = out_scen),
                 if (!inherits(obs_hat, "try-error")) {obs_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bart_hat, "try-error")) {bart_hat$tree_estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bart_hat, "try-error")) {bart_hat$smooth_estimate} else {rep(NA, length(a.vals))})
    
    #standard error
    se <- rbind(if (!inherits(obs_hat, "try-error")) {sqrt(obs_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(naive_hat, "try-error")) {sqrt(naive_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$tree_variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$smooth_variance)} else {rep(NA, length(a.vals))})
    
    # coverage probability
    cp <- rbind(if (!inherits(obs_hat, "try-error")) {as.numeric((est[2,] - 1.96*se[1,]) < est[1,] & (est[2,] + 1.96*se[1,]) > est[1,])} else {rep(NA, length(a.vals))},
                if (!inherits(naive_hat, "try-error")) {as.numeric((est[3,] - 1.96*se[2,]) < est[1,] & (est[3,] + 1.96*se[2,]) > est[1,])} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {as.numeric(bart_hat$hpdi[1,] < est[1,] & bart_hat$hpdi[2,] > est[1,])} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {as.numeric((est[5,] - 1.96*se[4,]) < est[1,] & (est[5,] + 1.96*se[4,]) > est[1,])} else {rep(NA, length(a.vals))})
    
    return(list(est = est, se = se, cp = cp))
     
  }, mc.cores = 25, mc.preschedule = TRUE)
  
  est <- abind(lapply(out, function(lst, ...) if (!inherits(lst, "try-error")) {lst$est} else {matrix(NA, ncol = length(a.vals), nrow = 6)}), along = 3)
  se <- abind(lapply(out, function(lst, ...) if (!inherits(lst, "try-error")) {lst$se} else {matrix(NA, ncol = length(a.vals), nrow = 5)}), along = 3)
  cp <- abind(lapply(out, function(lst, ...) if (!inherits(lst, "try-error")) {lst$cp} else {matrix(NA, ncol = length(a.vals), nrow = 5)}), along = 3)
  
  out_est <- t(apply(est, 1, rowMeans, na.rm = T))
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERF","DR","RC","BART","LOESS")
  
  compare <- matrix(est[1,,], nrow = length(a.vals), ncol = n.sim)
  out_bias <- t(apply(est[2:5,,], 1, function(x) rowMeans(abs(x - compare), na.rm = T)))
  colnames(out_bias) <- a.vals
  rownames(out_bias) <- c("DR","RC","BART","LOESS")
  
  out_sd <- t(apply(est[2:5,,], 1, function(x) apply(x, 1, sd, na.rm = T)))
  colnames(out_sd) <- a.vals
  rownames(out_sd) <- c("DR","RC","BART","LOESS")
  
  out_se <- t(apply(se, 1, rowMeans, na.rm = T))
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("DR","RC","BART","LOESS")
  
  out_cp <- t(apply(cp, 1, rowMeans, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("DR","RC","BART","LOESS")
  
  rslt <- list(scenario = scenario, est = out_est, bias = out_bias, sd = out_sd, se = out_se, cp = out_cp)
  filename <- paste0("~/Dropbox/Projects/ERC-EPE/Output/sim2/", paste(scenario, collapse = "_"),".RData")
  save(rslt, file = filename)
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.04)
sl.lib <- c("SL.mean","SL.glm")
n.sim <- 500

n <- c(400, 800)
mult <- c(5, 10)
sig_gps <- c(1, 2)
gps_scen <- c("a", "b")
out_scen <- c("a", "b")
pred_scen <- c("a", "b")

scen_mat <- expand.grid(n = n, mult = mult, sig_gps = sig_gps,  gps_scen = gps_scen, out_scen = out_scen, pred_scen = pred_scen, stringsAsFactors = FALSE)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- lapply(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib)
