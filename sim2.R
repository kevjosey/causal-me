### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(parallel)

# Code for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/gibbs-sampler.R")
source("D:/Github/causal-me/mclapply-hack.R")
source("D:/Github/causal-me/blp.R")
source("D:/Github/causal-me/erc.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # simulation arguments
  sig_gps <- 1
  sig_agg <- sqrt(2)
  sig_pred <- sqrt(0.5)
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  prob <- 0.1
  span <- 0.8
  
  # gen data arguments
  n <- scenario$n # c(500, 800)
  m <- n*scenario$mult# c(100, 200)
  
  # gibbs sampler stuff
  thin <- 20
  n.iter <- 1000
  n.adapt <- 100
  family <- poisson()
  
  # initialize output
  est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
  se <- array(NA, dim = c(n.sim, 4, length(a.vals)))
  
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
    
    stilde <- pred(s = s, star = star, w = w, sl.lib = sl.lib)
    a_w <- blp(s = stilde, s.id = s.id)
    a_x <- blp(s = stilde, s.id = s.id, x = x)
    
    w_new <- model.matrix(~ 0 + star*w)
    gibbs_w <- gibbs_dr(s = s, star = star, s.id = s.id, id = id, w = w_new,
                        n.iter = n.iter, n.adapt = n.adapt, thin = thin)
    
    gibbs_x <- gibbs_dr(s = s, star = star, s.id = s.id, id = id, w = w_new, x = x,
                        n.iter = n.iter, n.adapt = n.adapt, thin = thin)
    
    # blp w/ pred
    
    blp_w <- erc(y = y, a = a_w, x = x, offset = offset, family = family,
                 a.vals = a.vals, sl.lib = sl.lib, span = span)
    blp_x <- erc(y = y, a = a_x, x = x, offset = offset, family = family,
                 a.vals = a.vals, sl.lib = sl.lib, span = span)
    
    # gibbs w/ pred
    
    a_list_w <- split(gibbs_w$amat, seq(nrow(gibbs_w$amat)))
    a_list_x <- split(gibbs_x$amat, seq(nrow(gibbs_x$amat)))
    
    out_w <- mclapply.hack(1:length(a_list_x), function(k, ...){
      
      erc(y = y, a = a_list_w[[k]], x = x, offset = offset, family = family,
          a.vals = a.vals, sl.lib = sl.lib, span = span)
      
    }, mc.cores = 4)
    
    out_x <- mclapply.hack(1:length(a_list_x), function(k, ...){
      
      erc(y = y, a = a_list_x[[k]], x = x, offset = offset, family = family,
          a.vals = a.vals, sl.lib = sl.lib, span = span)
      
    }, mc.cores = 4)
    
    gibbs_est_w <- do.call(rbind, lapply(out_w, function(o) o$estimate))
    gibbs_var_w <- do.call(rbind, lapply(out_w, function(o) o$variance))
    
    gibbs_est_x <- do.call(rbind, lapply(out_x, function(o) o$estimate))
    gibbs_var_x <- do.call(rbind, lapply(out_x, function(o) o$variance))
    
    # estimates
    
    est[i,1,] <- predict_example(a = a.vals, x = x, out_scen = out_scen)
    est[i,2,] <- blp_w$estimate
    est[i,3,] <- colMeans(gibbs_est_w)
    est[i,4,] <- blp_x$estimate
    est[i,5,] <- colMeans(gibbs_est_x)
    
    #standard error
    
    var_w <- colMeans(gibbs_var_w) + (1 + 1/nrow(gibbs_est_w))*apply(gibbs_est_w, 2, var)
    var_x <- colMeans(gibbs_var_x) + (1 + 1/nrow(gibbs_est_x))*apply(gibbs_est_x, 2, var)
    
    se[i,1,] <- sqrt(blp_w$variance)
    se[i,2,] <- sqrt(var_w)
    se[i,3,] <- sqrt(blp_x$variance)
    se[i,4,] <- sqrt(var_x)
    
  }
  
  out_est <- colMeans(est, na.rm = T)
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("SAMPLE ERC","BLP W", "Gibbs W", "BLP X", "Gibbs X")
  
  out_se <- colMeans(se, na.rm = T)
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("SAMPLE ERC","BLP W", "Gibbs W", "BLP X", "Gibbs X")
  
  cp_blp_w <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,2,] - 1.96*se[i,1,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,2,] + 1.96*se[i,1,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  cp_gibbs_w <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,3,] - 1.96*se[i,2,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,3,] + 1.96*se[i,2,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  cp_blp_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,4,] - 1.96*se[i,3,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,4,] + 1.96*se[i,3,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  cp_gibbs_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,5,] - 1.96*se[i,4,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,5,] + 1.96*se[i,4,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  out_cp <- rbind(rowMeans(cp_blp_w, na.rm = T), rowMeans(cp_gibbs_w, na.rm = T),
                  rowMeans(cp_blp_x, na.rm = T), rowMeans(cp_gibbs_x, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("BLP W", "Gibbs W", "BLP X", "Gibbs X")
  
  return(list(est = out_est, se = out_se, cp = out_cp))
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.25)
sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction","SL.gam")
n.sim <- 1000

n <- c(500)
mult <- c(4)
gps_scen <- c("a", "b")
out_scen <- c("a", "b")

scen_mat <- expand.grid(n = n, mult = mult, gps_scen = gps_scen, out_scen = out_scen)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- mclapply.hack(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib, mc.cores = 3)
rslt <- list(est = est, scen_idx = scen_mat)

save(rslt, file = "D:/Github/causal-me/output/sim1_rslt.RData")