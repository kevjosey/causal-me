### Test DR estimator

rm(list = ls())

## Preliminaries

library(MASS)
library(cobalt)
library(SuperLearner)
library(data.table)

# Functions for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/hct-dr.R")

n.sim <- 100

# the true true
n <- 5000 # c(1000, 4000)
m <- 500 # c(100, 200)
sig_gps <- 2

a.vals <- seq(-1, 3, length.out = 21)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
span <- NULL
k <- 5
span.seq <- seq(0.15, 1, by = 0.05)

gps_scen <- "a"
out_scen <- "a"

est <- array(NA, dim = c(n.sim, 2, length(a.vals)))
sd <- matrix(NA, n.sim, length(a.vals))

for (i in 1:n.sim){
  
  print(i)
  
  dat <- gen_dr_data(n = n, m = m, sig_gps = sig_gps, gps_scen = gps_scen, out_scen = out_scen)
  
  y <- dat$y
  a <- dat$t 
  x <- dat$x
  y.id <- dat$id
  
  out <- hct_dr(y = y, a = a, x = x, y.id = y.id, a.vals = a.vals, span.seq = span.seq, k = k, sl.lib = sl.lib)
  est[i,1,] <- predict_example(a = a.vals, x = x, id = y.id, out_scen = out_scen)
  est[i,2,] <- out$estimate
  sd[i,] <- sqrt(out$variance)
  
}

cp <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*sd[i,]) < colMeans(est[,1,]) & (est[i,2,] + 1.96*sd[i,]) > colMeans(est[,1,])))

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
