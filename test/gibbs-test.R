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

# simulation arguments
n.sim <- 100
sig_epe <- sqrt(2)
sig_gps <- 1
out_scen <- "a"
gps_scen <- "a"

# gen data arguments
l <- 500 # c(500, 800)
m <- 100 # c(100, 200)
n <- 2000 # c(1000, 4000)

# dr arguments
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
span <- 0.9

# mcmc/prior values
shape <- 1e-5 # gamma shape
rate <- 1e-5 # gamma rate
scale <- 1e5 # normal scale
thin <- 10
n.adapt <- 1000
n.iter <- 10000

# initialize output
est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
se <- cp <- matrix(NA, n.sim, length(a.vals))

for (i in 1:n.sim){
  
  print(i)
  
  # generate data
  dat <- gen_data(l = l, m = m, n = n, sig_gps = sig_gps, sig_epe = sig_epe,
                  out_scen = out_scen, gps_scen = gps_scen)
  
  s.id <- dat$s.id
  y.id <- dat$y.id
  y <- dat$y
  s <- dat$s
  x <- dat$x
  a <- dat$a
  
  # remove clusters w/o exposure data
  ungroup <- which(!(1:m %in% unique(s.id)))
  if (length(ungroup)!= 0) {
    y <- dat$y[!(y.id %in% ungroup)]
    x <- dat$x[!(y.id %in% ungroup),]
    a_y <- dat$a_y[!(y.id %in% ungroup)]
    y.id <- dat$y.id[!(y.id %in% ungroup)]
  }
    
  fmla <- formula("~ x1 + x2 + x3 + x4")
  
  mcmc <- gibbs_dr(s = s, x = x, s.id = s.id, y.id = y.id, fmla = fmla,
                  shape = shape, rate = rate, scale = scale, 
                  thin = thin, n.iter = n.iter, n.adapt = n.adapt)
  
  bla <- blp(s = s, x = x, s.id = s.id, y.id = y.id, fmla = fmla)
  
  # multiple imputation
  a_list <- split(mcmc$amat_y, seq(nrow(mcmc$amat_y)))
  
  mi <- mclapply.hack(1:length(a_list), function(k, a_list, y, xmat, y.id, beta, sigma2, a.vals, sl.lib){
    
    hct_dr(y = y, a = a_list[[k]], x = xmat, y.id = y.id, a.vals = a.vals, k = 5,
           beta = beta[k,], sigma2 = sigma2[k], sl.lib = sl.lib)
    
  }, y = y, a_list = a_list, xmat = x, y.id = y.id, a.vals = a.vals,
  beta = mcmc$beta, sigma2 = mcmc$sigma2, sl.lib = sl.lib, mc.cores = 2)
  
  est_tmp <- matrix(unlist(lapply(mi, function(x) x$estimate)), ncol = length(mi))
  var_mat <- matrix(unlist(lapply(mi, function(x) x$variance)), ncol = length(mi))
  
  mi_est <- rowMeans(est_tmp)
  est_mat <- matrix(rep(mi_est, length(mi)), nrow = length(a.vals), ncol = length(mi))
  mi_var <- rowMeans(var_mat) + (1 + 1/length(mi))*rowSums((est_tmp - est_mat)^2)/(length(mi) - 1)
  
  t <- aggregate(s, by = list(s.id), mean)[,2]
  id <- unique(s.id)[order(unique(s.id))]
  t_y <- rep(NA, length(y.id))
  
  for (g in id) {
    t_y[y.id == g] <- t[id == g]
  }
  
  naive <- hct_dr(y = y, a = t_y, x = x, y.id = y.id, a.vals = a.vals, span = span, sl.lib = sl.lib)
  
  est[i,1,] <- predict_example(a = a.vals, x = x, id = y.id, out_scen = "a")
  est[i,2,] <- naive$estimate
  est[i,3,] <- out$est
  se[i,] <- out$se
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("TRUE SATE", "DR")

cp <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,]) > est[i,1,]))

out_cp <- rowMeans(cp, na.rm = T)
names(out_cp) <- a.vals

out_est
out_cp

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2,
      main = "Exposure Response Curve", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "blue", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "green", lwd = 2)
legend(2, 0.2, legend=c("Sample ERC", "Naive", "SIMEX"), col=c("red", "blue", "green"), lwd=2, cex=0.8)
