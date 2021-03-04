### Test DR estimator

rm(list = ls())

## Preliminaries

library(MASS)
library(cobalt)
library(SuperLearner)
library(data.table)
library(parallel)

# Code for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/simex-dr.R")
source("D:/Github/causal-me/hct-dr.R")
source("D:/Github/causal-me/mclapply-hack.R")

n.sim <- 100
n.boot <- 200

# the true true
# the true true
m <- 500 # c(100, 200)
n <- 5000 # c(1000, 4000)

# True values
sig_epe <- 3
sig_gps <- 5
prob <- 0.1

# other preliminaries
span.seq <- seq(0.5, 1, by = 0.05)
span <- 0.8
k <- 5
degree <- 2
lambda <- seq(0.1, 2.1, by = 0.25)

a.vals <- seq(-1, 3, by = 0.5)

est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
se <- cp <- matrix(NA, n.sim, length(a.vals))

for (i in 1:n.sim){
  
  print(i)
  
  dat <- gen_simex_data(m = m, n = n, sig_epe = sig_epe, sig_gps = sig_gps)
  
  id <- dat$id
  y.id <- dat$y.id
  
  y <- dat$y
  a <- dat$a
  z <- dat$z
  x <- dat$x
  z_y <- dat$z_y
  a_y <- dat$a_y
  
  out <- simex_dr(z = z, y = y, x = x, id = id, y.id = y.id, sigma = sig_epe, a.vals = a.vals,
                  n.boot = n.boot, degree = degree, lambda = lambda, span = span, span.seq = span.seq, k = k)
  
  naive <- hct_dr(y = y, a = z_y, x = x, y.id = y.id, a.vals = a.vals, span.seq = span.seq, k = k, sl.lib = sl.lib)
  
  est[i,1,] <- predict_example(a.vals = a.vals, x = x, y.id = y.id, out_scen = "a")
  est[i,2,] <- naive$estimate
  est[i,3,] <- out$estimate
  se[i,] <- sqrt(out$variance)
  
}

cp <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,]) > est[i,1,]))

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("TRUE SATE", "NAIVE", "SIMEX")
out_cp <- rowMeans(cp, na.rm = T)
names(out_cp) <- a.vals

out_est
out_cp

plot(a.vals, est[i,3,], type = "l", col = "grey", lwd = 2, main = "Exposure Response Curve", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))

for(i in 2:nrow(cp))
  lines(a.vals, est[i,3,], type = "l", col = "grey", lwd = 2)

lines(a.vals, colMeans(est[,1,], na.rm = T), type = "l", col = "red", lwd = 2)
lines(a.vals, colMeans(est[,2,], na.rm = T), type = "l", col = "blue", lwd = 2)
