### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(earth)
library(parallel)

# Code for generating and fitting data
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/gen-data.R")
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/gibbs-sampler.R")
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/blp.R")
source("~/shared_space/ci3_analysis/josey_causal_me/causal-me/erc.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # simulation arguments
  sig_gps <- 1
  sig_agg <- scenario$sig_agg
  sig_pred <- scenario$sig_pred
  prob <- scenario$prob
  span <- 0.75
  gps_scen <- "a"
  out_scen <- "a"
  
  # gen data arguments
  n <- scenario$n
  m <- scenario$mult*n
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
    
    # exposure predictions
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
                           a.vals = a.vals, sl.lib = sl.lib, span = span), silent = TRUE)
    naive_hat <- try(erc(y = y, a = z_hat, x = x, offset = offset, family = family,
                         a.vals = a.vals, sl.lib = sl.lib, span = span), silent = TRUE)
    
    # blp approach
    blp_tilde <- try(erc(y = y, a = a_tilde, x = x, offset = offset, family = family,
                         a.vals = a.vals, sl.lib = sl.lib, span = span), silent = TRUE)
    blp_hat <- try(erc(y = y, a = a_hat, x = x, offset = offset, family = family,
                       a.vals = a.vals, sl.lib = sl.lib, span = span), silent = TRUE)
    
    # estimates
    est[i,1,] <- predict_example(a = a.vals, x = x, out_scen = out_scen)
    est[i,2,] <- if (!inherits(naive_tilde, "try-error")) {naive_tilde$estimate} else {rep(NA, length(a.vals))}
    est[i,3,] <- if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))}
    est[i,4,] <- if (!inherits(blp_tilde, "try-error")) {blp_tilde$estimate} else {rep(NA, length(a.vals))}
    est[i,5,] <- if (!inherits(blp_hat, "try-error")) {blp_hat$estimate} else {rep(NA, length(a.vals))}
    
    # standard errors
    se[i,1,] <- if (!inherits(naive_tilde, "try-error")) {sqrt(naive_tilde$variance)} else {rep(NA, length(a.vals))}
    se[i,2,] <- if (!inherits(naive_hat, "try-error")) {sqrt(naive_hat$variance)} else {rep(NA, length(a.vals))}
    se[i,3,] <- if (!inherits(blp_tilde, "try-error")) {sqrt(blp_tilde$variance)} else {rep(NA, length(a.vals))}
    se[i,4,] <- if (!inherits(blp_hat, "try-error")) {sqrt(blp_hat$variance)} else {rep(NA, length(a.vals))}
    
  }
  
  out_est <- colMeans(est, na.rm = T)
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERC", "Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  out_se <- colMeans(se, na.rm = T)
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  cp_naive_w <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,2,] - 1.96*se[i,1,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,2,] + 1.96*se[i,1,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  cp_naive_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,3,] - 1.96*se[i,2,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,3,] + 1.96*se[i,2,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  cp_blp_w <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,4,] - 1.96*se[i,3,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,4,] + 1.96*se[i,3,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  cp_blp_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[i,5,] - 1.96*se[i,4,]) < colMeans(est[,1,],na.rm = TRUE) & 
                 (est[i,5,] + 1.96*se[i,4,]) > colMeans(est[,1,],na.rm = TRUE)))
  
  out_cp <- rbind(rowMeans(cp_naive_w, na.rm = T), rowMeans(cp_naive_x, na.rm = T),
                  rowMeans(cp_blp_w, na.rm = T), rowMeans(cp_blp_x, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("Naive Tilde", "Naive Hat", "BLP Tilde", "BLP Hat")
  
  return(list(est = out_est, se = out_se, cp = out_cp))
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.1)
sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction","SL.earth")
n.sim <- 1000

n <- c(200, 500)
mult <- c(5, 10)
sig_pred <- c(0, sqrt(0.5)) 
sig_agg <- c(0, sqrt(2))
prob <- c(0.1, 0.2)

scen_mat <- expand.grid(n = n, mult = mult, sig_agg = sig_agg, sig_pred = sig_pred, prob = prob)
scen_mat <- round(scen_mat, 3)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- mclapply(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib, mc.cores = 16)
rslt <- list(est = est, scen_idx = scen_mat)

save(rslt, file = "~/shared_space/ci3_analysis/josey_causal_me/Output/sim1_rslt.RData")

for (k in 1:length(rslt$est)){
  
  filename <- paste0("~/shared_space/ci3_analysis/josey_causal_me/Output/", paste(rslt$scen_idx[k,], collapse = "_"), ".pdf")
  pdf(file = filename)
  plotname <- paste(rslt$scen_idx[k,], collapse = "_")
  
  plot(a.vals, rslt$est[[k]]$est[1,], type = "l", col = "darkgreen", lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event",
       ylim = c(0,0.1))
  lines(a.vals, rslt$est[[k]]$est[2,], type = "l", col = "red", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[k]]$est[3,], type = "l", col = "blue", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[k]]$est[4,], type = "l", col = "red", lwd = 2, lty = 3)
  lines(a.vals, rslt$est[[k]]$est[5,], type = "l", col = "blue", lwd = 2, lty = 3)
  
  legend(6, 0.1, legend=c("True ERC", "Without Prediction Correction",
                          "With Prediction Correction", "Without Aggregation Correction",
                          "With Aggregation Correction"),
         col=c("darkgreen", "red", "blue", "black", "black"),
         lty = c(1,1,1,2,3), lwd=2, cex=0.8)
  
  dev.off()
  
}
