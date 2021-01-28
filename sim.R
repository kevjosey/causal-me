### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(parallel)

# Code for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/gibbs-dr.R")
source("D:/Github/causal-me/simex-dr.R")
source("D:/Github/causal-me/hct-dr.R")
source("D:/Github/causal-me/mclapply-hack.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # gen_data arguments
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  l <- scenario$l
  m <- scenario$m
  n <- scenario$n
  
  # for now these are fixed
  sig_epe <- sqrt(2)
  sig_gps <- 1
  
  # mcmc/prior arguments
  shape <- 1e-5 # gamma shape
  rate <- 1e-5 # gamma rate
  scale <- 1e5 # normal scale
  thin <- 20
  n.adapt <- 1000
  n.iter <- 1000
  
  # simex arguments
  n.boot <- 50
  degree <- 2
  lambda <- seq(0.1, 2.1, by = 0.2)
  
  # initialize output
  est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
  se <- array(NA, dim = c(n.sim, 4, length(a.vals)))

  for (i in 1:n.sim){
    
    print(i)
    
    # generate data
    dat <- gen_data(l = l, m = m, n = n, 
                    sig_gps = sig_gps, sig_epe = sig_epe, 
                    gps_scen = gps_scen, out_scen = out_scen)
    
    s.id <- dat$s.id
    y.id <- dat$y.id
    y <- dat$y
    s <- dat$s
    x <- dat$x
    a <- dat$a
    
    # remove any clusters w/o exposure data (random process)
    ungroup <- which(!(1:m %in% unique(s.id)))
    if (length(ungroup)!= 0) {
      y <- dat$y[!(y.id %in% ungroup)]
      x <- dat$x[!(y.id %in% ungroup),]
      y.id <- dat$y.id[!(y.id %in% ungroup)]
    }
    
    id <- unique(y.id)[order(unique(y.id))]
    
    fmla.s <- formula("~ w1 + w2 + w3 + w4")
    fmla.a <- formula("~ x1 + x2 + x3 + x4")
    
    # estimate parameters and generate latent variable
    mcmc <- gibbs_dr(s = s, y = y, x = x, s.id = s.id, y.id = y.id, 
                    fmla.a = fmla.a, n.iter = n.iter, n.adapt = n.adapt,
                    shape = shape, rate = rate, scale = scale, thin = thin)
    
    z_list <- split(mcmc$zMat_y, seq(nrow(mcmc$zMat_y)))
    
    # multiple imputation
    mi <- mclapply.hack(z_list, function(z, y, xmat, y.id, a.vals, sl.lib){
  
      hct_dr(y = y, a = z, x = xmat, y.id = y.id, a.vals = a.vals, k = 5, sl.lib = sl.lib)
  
    }, y = y, xmat = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, mc.cores = 2)
    
    est_tmp <- matrix(unlist(lapply(mi, function(x) x$estimate)), ncol = length(mi))
    var_mat <- matrix(unlist(lapply(mi, function(x) x$variance)), ncol = length(mi))
    
    mi_est <- rowMeans(est_tmp)
    est_mat <- matrix(rep(mi_est, length(mi)), nrow = length(a.vals), ncol = length(mi))
    mi_var <- rowMeans(var_mat) + (1 + 1/length(mi))*rowSums((est_tmp - est_mat)^2)/(length(mi) - 1)
    
    # simulation extrapolation
    simex <- simex_dr(z = colMeans(mcmc$zMat), y = y, x = x, id = id, y.id = y.id, 
                      sigma = sqrt(mean(mcmc$tau2)/table(s.id)), n.boot = n.boot, degree = degree,
                      lambda = lambda, a.vals = a.vals, k = 5, sl.lib = sl.lib, mc.cores = 2)
    
    simex_est <- simex$estimate
    simex_var <- simex$variance
    
    # single imputation
    si <- hct_dr(y = y, a = colMeans(mcmc$zMat_y), x = x, y.id = y.id, a.vals = a.vals, k = 5, sl.lib = sl.lib)
    
    si_est <- si$estimate
    si_var <- si$variance
    
    # naive approach
    t <- aggregate(s, by = list(s.id), mean)[,2]
    id <- unique(s.id)[order(unique(s.id))]
    a_y <- rep(NA, length(y.id))
    
    for (g in id)
      a_y[y.id == g] <- t[id == g]
    
    naive <- hct_dr(y = y, a = a_y, x = x, y.id = y.id, a.vals = a.vals, k = 5, sl.lib = sl.lib)
    
    naive_est <- naive$estimate
    naive_var <- naive$variance
    
    # combine output
    est[i,1,] <- predict_example(a.vals = a.vals, x = x, y.id = y.id, out_scen = out_scen)
    est[i,2:5,] <- rbind(naive_est, si_est, simex_est, mi_est)
    se[i,1:4,] <- sqrt(rbind(naive_var, si_var, simex_var, mi_var))
    
  }

  out_est <- colMeans(est, na.rm = T)
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("TRUE SATE", "NAIVE", "SI", "SIMEX", "MI")
  
  out_se <- colMeans(se, na.rm = T)
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("NAIVE", "SI", "SIMEX", "MI")
  
  cp_list <- lapply(1:n.sim, function(i,...)
    list(naive = as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]),
          si = as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]),
          simex = as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]),
          mi = as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & (est[i,5,] + 1.96*se[i,4,]) > est[i,1,])))
  
  cp_naive <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$naive)), ncol = length(cp_list)), na.rm = TRUE)
  cp_si <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$si)), ncol = length(cp_list)), na.rm = TRUE)
  cp_simex <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$simex)), ncol = length(cp_list)), na.rm = TRUE)
  cp_mi <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$mi)), ncol = length(cp_list)), na.rm = TRUE)
  
  out_cp <- rbind(NAIVE = cp_naive, SI = cp_si, SIMEX = cp_simex, MI = cp_mi)
  colnames(out_cp) <- a.vals

  out <- list(est = out_est, se = out_se, cp = out_cp)
  return(out)
  
}
  
# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
n.sim <- 100

n <- c(2000)
m <- c(100)
l <- c(500)
out_scen <- c("a", "b")
gps_scen <- c("a", "b")

scen_mat <- expand.grid(n = n, m = m, l = l, out_scen = out_scen, gps_scen = gps_scen)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
rslt <- mclapply.hack(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib, mc.cores = 2)
rslt <- list(unlist(rslt), scen_mat)

save(rslt, file = "D:/Github/causal-me/output/rslt.RData")

pdf(file = "D:/Github/causal-me/output/ERC-plot.pdf")
par(mfrow = c(2, 2))
for (k in 1:length(rslt)) {
  
  plot(a.vals, rslt[[k]]$est[1,], type = "l", col = "red", lwd = 2,
       main = paste("scenario", k), xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
  lines(a.vals, rslt[[k]]$est[2,], type = "l", col = "orange", lwd = 2)
  lines(a.vals, rslt[[k]]$est[3,], type = "l", col = "green", lwd = 2)
  lines(a.vals, rslt[[k]]$est[4,], type = "l", col = "blue", lwd = 2)
  lines(a.vals, rslt[[k]]$est[5,], type = "l", col = "purple", lwd = 2)

  if(k == length(rslt))
    legend(2, 0.2, legend=c("Sample ERC", "Naive", "SI", "SIMEX", "MI"), col=c("red", "orange", "green", "blue", "purple"), lwd=2, cex=0.8)
  
}
dev.off()
