### Test DR estimator

rm(list = ls())

## Preliminaries

library(MASS)
library(cobalt)
library(SuperLearner)
library(data.table)
library(parallel)

# Code for generating and fitting data
source("D:/Dropbox (Personal)/ERC-EPE/Code/gen-data.R")
source("D:/Dropbox (Personal)/ERC-EPE/Code/simex-dr.R")
source("D:/Dropbox (Personal)/ERC-EPE/Code/hct-dr.R")

n.sim <- 100
n.boot <- 20

# the true true
# the true true
m <- 100 # c(100, 200)
n <- 2000 # c(1000, 4000)

sig_epe <- 2
sig_gps <- 1
prob <- 0.1

# ml libraries
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
span.seq <- seq(0.5, 1, by = 0.05)
span <- 0.8
k <- 5
degree <- 2
lambda <- seq(0.1, 2.1, by = 0.1)

a.vals <- seq(-1, 3, by = 0.5)

est <- array(NA, dim = c(n.sim, 2, length(a.vals)))
se <- cp <- matrix(NA, n.sim, length(a.vals))

for (i in 1:n.sim){
  
  print(i)
  
  dat <- gen_agg_data(m = m, n = n, prob = prob)
  
  s.id <- dat$s.id
  y.id <- dat$y.id
  
  ungroup <- which(!(1:m %in% unique(s.id)))
  
  y <- dat$y
  a <- dat$a
  z <- dat$z
  x <- dat$x
  
  if (length(ungroup)!= 0) {
    y <- dat$y[!(y.id %in% ungroup)]
    x <- dat$x[!(y.id %in% ungroup),]
    y.id <- dat$y.id[!(y.id %in% ungroup)]
  }
  
  out <- simex_dr(z = z, y = y, x = x, s.id = s.id, y.id = y.id, sig_epe = sig_epe, n.boot = n.boot, degree = degree,
                  a.vals = a.vals, lambda = lambda, span = span, span.seq = span.seq, k = k, sl.lib = sl.lib)
  
  est[i,1,] <- predict.example(a = a.vals, x = x, id = y.id)
  est[i,2,] <- out$estimate
  se[i,] <- sqrt(out$variance)
  
}

cp <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,]) > est[i,1,]))

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("TRUE SATE", "DR")
out_cp <- rowMeans(cp, na.rm = T)
names(out_cp) <- a.vals

out_est
out_cp

plot(a.vals, est[i,2,], type = "l", col = "grey", lwd = 2, main = "Exposure Response Curve", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))

for(i in 2:nrow(cp))
  lines(a.vals, est[i,2,], type = "l", col = "grey", lwd = 2)

lines(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2)
