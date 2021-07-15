### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(earth)
library(parallel)
library(abind)

# Code for generating and fitting data
source("~/Github/causal-me/mclapply-hack.R")
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/erc.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # simulation arguments
  sig_gps <- 1
  sig_agg <- scenario$sig_agg
  sig_pred <- scenario$sig_pred
  gps_scen <- "a"
  out_scen <- "a"
  pred_scen <- "b"
  
  # gen data arguments
  n <- scenario$n
  mult <- scenario$mult
  prob <- scenario$prob
  
  # dr arguments
  family <- poisson()
  span <- ifelse(n == 800, 0.125, 0.25)
  deg.num <- 2

  # initialize output
  est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
  se <- array(NA, dim = c(n.sim, 4, length(a.vals)))
  
  print(scenario)
  
  # generate data
  dat_list <- replicate(n.sim, gen_data(n = n, mult = mult, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
                  pred_scen = pred_scen, out_scen = out_scen, gps_scen = gps_scen))
  
  out <- mclapply.hack(1:n.sim, function(i, ...){
    
    dat <- dat_list[,i]
    
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
    
    s_hat <- pred(s = s, star = s_tilde, w = w, sl.lib = sl.lib)
    a_tilde <- blp(s = s_tilde, s.id = s.id, x = x)
    a_hat <- blp(s = s_hat, s.id = s.id, x = x)
    
    z_tilde_tmp <- aggregate(s_tilde, by = list(s.id), mean)
    z_hat_tmp <- aggregate(s_hat, by = list(s.id), mean)
    
    z_hat <- z_tilde <- rep(NA, length(id))
    
    for (g in id) {
      
      z_tilde[id == g] <- z_tilde_tmp[z_tilde_tmp[,1] == g,2]
      z_hat[id == g] <- z_hat_tmp[z_hat_tmp[,1] == g,2]
      
    }
    
    # naive
    naive_tilde <- try(erc(y = y, a = z_tilde, x = x, offset = offset, family = family,
                           a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    naive_hat <- try(erc(y = y, a = z_hat, x = x, offset = offset, family = family,
                         a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    
    # blp approach
    blp_tilde <- try(erc(y = y, a = a_tilde, x = x, offset = offset, family = family,
                         a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    blp_hat <- try(erc(y = y, a = a_hat, x = x, offset = offset, family = family,
                       a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    
    # estimates
    est <- rbind(predict_example(a = a.vals, x = x, out_scen = out_scen),
                 if (!inherits(naive_tilde, "try-error")) {naive_tilde$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(blp_tilde, "try-error")) {blp_tilde$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(blp_hat, "try-error")) {blp_hat$estimate} else {rep(NA, length(a.vals))})
    
    # standard errors
    se <- rbind(if (!inherits(naive_tilde, "try-error")) {sqrt(naive_tilde$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(naive_hat, "try-error")) {sqrt(naive_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(blp_tilde, "try-error")) {sqrt(blp_tilde$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(blp_hat, "try-error")) {sqrt(blp_hat$variance)} else {rep(NA, length(a.vals))})
    
    cp <- rbind(if (!inherits(naive_tilde, "try-error")) {as.numeric((est[2,] - 1.96*se[1,]) < est[1,] & (est[2,] + 1.96*se[1,]) > est[1,])} 
                else {rep(NA, length(a.vals))},
                if (!inherits(naive_hat, "try-error")) {as.numeric((est[3,] - 1.96*se[2,]) < est[1,] & (est[3,] + 1.96*se[2,]) > est[1,])}
                else {rep(NA, length(a.vals))},
                if (!inherits(blp_tilde, "try-error")) {as.numeric((est[4,] - 1.96*se[3,]) < est[1,] & (est[4,] + 1.96*se[3,]) > est[1,])}
                else {rep(NA, length(a.vals))},
                if (!inherits(blp_hat, "try-error")) {as.numeric((est[5,] - 1.96*se[4,]) < est[1,] & (est[5,] + 1.96*se[4,]) > est[1,])}
                else {rep(NA, length(a.vals))})
    
    return(list(est = est, se = se, cp = cp))
    
  }, mc.cores = 8)
  
  est <- abind(lapply(out, function(lst, ...) lst$est), along = 3)
  se <- abind(lapply(out, function(lst, ...) lst$se), along = 3)
  cp <- abind(lapply(out, function(lst, ...) lst$cp), along = 3)
  
  out_est <- t(apply(est, 1, rowMeans, na.rm = T))
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERF", "Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  compare <- matrix(est[1,,], nrow = length(a.vals), ncol = n.sim)
  out_bias <- t(apply(est[2:5,,], 1, function(x) rowMeans(abs(x - compare), na.rm = T)))
  colnames(out_bias) <- a.vals
  rownames(out_bias) <- c("Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  out_sd <- t(apply(est[2:5,,], 1, function(x) apply(x, 1, sd, na.rm = T)))
  colnames(out_sd) <- a.vals
  rownames(out_sd) <- c("Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  out_se <- t(apply(se, 1, rowMeans, na.rm = T))
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  out_cp <- t(apply(cp, 1, rowMeans, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  return(list(est = out_est, bias = out_bias, sd = out_sd, se = out_se, cp = out_cp))
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.04)
sl.lib <- c("SL.mean","SL.glm")
n.sim <- 1000

n <- c(400, 800)
mult <- c(5, 10)
sig_pred <- c(0, sqrt(0.5)) 
sig_agg <- c(0, sqrt(2))
prob <- c(0.1, 0.2)

scen_mat <- expand.grid(n = n, mult = mult, sig_agg = sig_agg, sig_pred = sig_pred, prob = prob)
scen_mat <- round(scen_mat, 3)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- lapply(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib)
rslt <- list(est = est, scen_idx = scen_mat)

save(rslt, file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/sim_1.RData")