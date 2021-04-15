### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(splines)
library(parallel)

# Code for generating and fitting data
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/gen-data.R")
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/gibbs-sampler.R")
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/blp.R")
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/erc.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # simulation arguments
  sig_gps <- 1
  sig_agg <- sqrt(2)
  sig_pred <- sqrt(0.5)
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  prob <- 0.2
  span <- 0.75
  
  # gen data arguments
  n <- scenario$n # c(500, 800)
  m <- n*scenario$mult # c(100, 200)
  
  # gibbs sampler stuff
  thin <- 50
  n.iter <- 5000
  n.adapt <- 500
  family <- poisson()
  
  # initialize output
  est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
  se <- array(NA, dim = c(n.sim, 2, length(a.vals)))
  
  for (i in 1:n.sim){
    
    print(i)
    
    # generate data
    dat <- gen_data(m = m, n = n, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
                    pred_scen = "b", out_scen = out_scen, gps_scen = gps_scen)
    
    # zipcode index
    s.id <- dat$s.id
    id <- dat$id
    offset <- log(dat$offset)
    
    # data
    y <- dat$y
    x <- dat$x
    a <- dat$a
    w <- dat$w
    star <- dat$star
    
    # validation subset
    s <- dat$s*rbinom(m, 1, 0.1)
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
    s_hat <- pred(s = s, star = star, w = w, sl.lib = sl.lib)
    a_x <- blp(s = s_hat, s.id = s.id, x = x)
    
    # blp w/ pred
    blp_x <- try(erc(y = y, a = a_x, x = x, offset = offset, family = family,
                     a.vals = a.vals, sl.lib = sl.lib, span = span), silent = TRUE)
    
    # Bayesian Approach
    gibbs_x <- try(gibbs_dr(s = s, star = star, y = y, offset = offset,
                            s.id = s.id, id = id, w = w, x = x, family = family,
                            n.iter = n.iter, n.adapt = n.adapt, thin = thin, 
                            h.a = 1, h.gamma = 0.3, deg.num = 4,
                            a.vals = a.vals, span = span, mc.cores = 1), silent = TRUE)
    
    # estimates
    est[i,1,] <- predict_example(a = a.vals, x = x, out_scen = out_scen)
    est[i,2,] <- if (!inherits(blp_x, "try-error")) {blp_x$estimate} else {rep(NA, length(a.vals))}
    est[i,3,] <- if (!inherits(gibbs_x, "try-error")) {gibbs_x$estimate} else {rep(NA, length(a.vals))}
    
    #standard error
    se[i,1,] <- if (!inherits(blp_x, "try-error")) {sqrt(blp_x$variance)} else {rep(NA, length(a.vals))}
    se[i,2,] <- if (!inherits(gibbs_x, "try-error")) {sqrt(gibbs_x$variance)} else {rep(NA, length(a.vals))}
    
  }
  
  out_est <- colMeans(est, na.rm = T)
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERC","BLP","Bayes")
  
  out_se <- colMeans(se, na.rm = T)
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("ERC","BLP","Bayes")
  
  cp_blp_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,2,] - 1.96*se[i,1,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,2,] + 1.96*se[i,1,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  cp_gibbs_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,3,] - 1.96*se[i,2,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,3,] + 1.96*se[i,2,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  out_cp <- rbind(rowMeans(cp_blp_x, na.rm = T), rowMeans(cp_gibbs_x, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("BLP", "Bayes")
  
  return(list(est = out_est, se = out_se, cp = out_cp))
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.1)
sl.lib <- c("SL.glm")
n.sim <- 200

n <- c(200,500)
mult <- c(5,10)
gps_scen <- c("a", "b")
out_scen <- c("a", "b")

scen_mat <- expand.grid(n = n, mult = mult, gps_scen = gps_scen, out_scen = out_scen)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- mclapply(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib, mc.cores = 16)
rslt <- list(est = est, scen_idx = scen_mat)

save(rslt, file = "~/shared_space/ci3_analysis/josey_causal_me/Output/sim2_rslt.RData")