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
prob <- 0.1

# gen data arguments
l <- 1000 # c(500, 800)
m <- 200 # c(100, 200)
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
n.adapt <- 100
n.iter <- 1000

# simex arguments
n.boot <- 100
degree <- 2
lambda <- seq(0.1, 2.1, by = 0.25)

# initialize output
est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
se <- array(NA, dim = c(n.sim, 4, length(a.vals)))

for (i in 1:n.sim){
  
  print(i)
  
  # generate data
  dat <- gen_bayes_data(l = l, m = m, n = n, sig_gps = sig_gps, sig_epe = sig_epe, prob = prob)
  
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
    y.id <- dat$y.id[!(y.id %in% ungroup)]
  }
  
  id <- unique(y.id)[order(unique(y.id))]
  
  fmla.s <- formula("~ w1 + w2 + w3 + w4")
  fmla.a <- formula("~ x1 + x2 + x3 + x4")
  
  mcmc <- gibbs_dr(s = s, y = y, x = x, s.id = s.id, y.id = y.id, 
                  fmla.s = fmla.s, fmla.a = fmla.a, n.iter = n.iter, n.adapt = n.adapt,
                  shape = shape, rate = rate, scale = scale, thin = thin)
  
  # multiple imputation
  mi <- mclapply.hack(1:nrow(mcmc$zMat), function(k, ...){

    hct_dr(y = y, a = mcmc$zMat_y[k,], x = x, y.id = y.id, a.vals = a.vals, span = span, sl.lib = sl.lib)

  })

  mi_est <- colMeans(t(matrix(unlist(lapply(mi, function(x) x$estimate)), ncol = length(mi))))
  mi_var <- colMeans(t(matrix(unlist(lapply(mi, function(x) x$variance)), ncol = length(mi))))
  
  # simulation extrapolation

  simex <- simex_dr(z = colMeans(mcmc$zMat), y = y, x = x, id = id, y.id = y.id, 
                    sigma = sqrt(mean(sigma2)/table(s.id)), n.boot = n.boot, degree = degree,
                    lambda = lambda, a.vals = a.vals, span = span, span.seq = span.seq, k = k, sl.lib = sl.lib)
  
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
  
  naive <- hct_dr(y = y, a = a_y, x = x, y.id = y.id, a.vals = a.vals, span = span, sl.lib = sl.lib)
  
  naive_est <- naive$estimate
  naive_var <- naive$variance
  
  # combine output
  est[i,1,] <- predict_example(a = a.vals, x = x, id = y.id, out_scen = "a")
  est[i,2:5,] <- rbind(naive_est, si_est, simex_est, mi_est)
  se[i,1:4,] <- sqrt(rbind(naive_var, si_var, simex_var, mi_var))
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("TRUE SATE", "NAIVE", "SI", "SIMEX", "MI")

cp <- sapply(1:n.sim, function(i,...)
  cbind(naive = as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]),
        naive = as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]),
        naive = as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]),
        naive = as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & (est[i,5,] + 1.96*se[i,4,]) > est[i,1,])))

out_cp <- rowMeans(cp, na.rm = T)
colnames(out_cp) <- a.vals

out_est
out_cp

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2,
     main = "Exposure Response Curve", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "blue", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "green", lwd = 2)
legend(2, 0.2, legend=c("Sample ERC", "Naive", "SIMEX"), col=c("red", "blue", "green"), lwd=2, cex=0.8)
