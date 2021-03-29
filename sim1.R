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
  sig_agg <- scenario$sig_agg
  sig_pred <- scenario$sig_pred
  prob <- scenario$prob
  span <- 0.8
  gps_scen <- "a"
  out_scen <- "a"
  
  # gen data arguments
  n <- scenario$n
  m <- scenario$mult*n
  family <- poisson()
  
  # initialize output
  est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
  
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
    s_tilde <- dat$star
    
    # validation subset
    s <- dat$s*rbinom(m, 1, prob)
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
    
    a_tilde <- blp(s = s_tilde, s.id = s.id)
    a_hat <- blp(s = s_hat, s.id = s.id)
    
    z_tilde_tmp <- aggregate(s_tilde, by = list(s.id), mean)
    z_hat_tmp <- aggregate(s_hat, by = list(s.id), mean)
    
    z_hat <- z_tilde <- rep(NA, length(id))
    
    for (g in id) {
      
      z_tilde[id == g] <- z_tilde_tmp[z_tilde_tmp[,1] == g,2]
      z_hat[id == g] <- z_hat_tmp[z_hat_tmp[,1] == g,2]
      
    }
    
    # naive
    
    naive_tilde <- erc(y = y, a = z_tilde, x = x, offset = offset, family = family,
                       a.vals = a.vals, sl.lib = sl.lib, span = span)
    naive_hat <- erc(y = y, a = z_hat, x = x, offset = offset, family = family,
                     a.vals = a.vals, sl.lib = sl.lib, span = span)
    
    # blp w/ pred
    
    blp_tilde <- erc(y = y, a = a_tilde, x = x, offset = offset, family = family,
                     a.vals = a.vals, sl.lib = sl.lib, span = span)
    blp_hat <- erc(y = y, a = a_hat, x = x, offset = offset, family = family,
                   a.vals = a.vals, sl.lib = sl.lib, span = span)
    
    # estimates
    
    est[i,1,] <- predict_example(a = a.vals, x = x, out_scen = out_scen)
    est[i,2,] <- naive_tilde$estimate
    est[i,3,] <- naive_hat$estimate
    est[i,4,] <- blp_tilde$estimate
    est[i,5,] <- blp_hat$estimate
    
  }
  
  out_est <- colMeans(est, na.rm = T)
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("SAMPLE ERC","Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")

  return(out_est)
  
}
  
# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.25)
sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction")
n.sim <- 1000

n <- c(200, 500)
mult <- c(2, 4)
sig_pred <- c(sqrt(0.5), sqrt(2)) 
sig_agg <- c(sqrt(0.5), sqrt(2))
prob <- c(0.1, 0.2)

scen_mat <- expand.grid(n = n, mult = mult, sig_agg = sig_agg, sig_pred = sig_pred, prob = prob)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- mclapply.hack(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib, mc.cores = 3)
rslt <- list(est = est, scen_idx = scen_mat)

save(rslt, file = "D:/Github/causal-me/output/sim1_rslt.RData")

for (k in 1:length(rslt$est)){

  filename <- paste0("D:/Github/causal-me/output/ERC_", paste(rslt$scen_idx[k,], collapse = "_"), ".pdf")
  pdf(file = filename)
  plotname <- paste(rslt$scen_idx[k,], collapse = "_")

  plot(a.vals, rslt$est[[k]][1,], type = "l", col = "darkgreen", lwd = 2,
       main = "Exposure = a, Outcome = a", xlab = "Exposure", ylab = "Rate of Event",
       ylim = c(0,0.1))
  lines(a.vals, rslt$est[[k]][2,], type = "l", col = "red", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[k]][3,], type = "l", col = "blue", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[k]][4,], type = "l", col = "red", lwd = 2, lty = 3)
  lines(a.vals, rslt$est[[k]][5,], type = "l", col = "blue", lwd = 2, lty = 3)

  legend(6, 0.1, legend=c("Sample ERC", "Without Prediction Correction",
                          "With Prediction Correction", "Without Aggregation Correction",
                          "With Aggregation Correction"),
         col=c("darkgreen", "red", "blue", "black", "black"),
         lty = c(1,1,1,2,3), lwd=2, cex=0.8)

  dev.off()
  
}

