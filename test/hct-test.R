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
n <- 1000 # c(1000, 4000)
m <- 100 # c(100, 200)
sig_gps <- 2

a.vals <- seq(-1, 3, length.out = 21)
span <- NULL
k <- 5
span.seq <- seq(0.15, 1, by = 0.05)

gps_scen <- "a"
out_scen <- "b"

est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
se <- array(NA, dim = c(n.sim, 2, length(a.vals)))

for (i in 1:n.sim){
  
  print(i)
  
  offset <- NULL
  
  dat <- gen_dr_data(n = n, m = m, sig_gps = sig_gps, gps_scen = gps_scen, out_scen = out_scen)
  
  y <- dat$y
  a <- dat$a
  a_y <- dat$a_y
  x <- dat$x
  y.id <- dat$id
  
  out1 <- hct_dr(y = y, a = a_y, x = x, a.vals = a.vals, y.id = y.id, span.seq = span.seq, k = k)
  
  est[i,1,] <- predict_example(a.vals = a.vals, x = x, y.id = y.id, out_scen = out_scen)
  est[i,2,] <- out1$estimate
  se[i,1,] <- sqrt(out1$variance)
  
  df.tmp <- setDT(data.frame(y.id = y.id, y = y, a = a_y, x = x))
  df <- df.tmp[,lapply(.SD, mean), by = y.id][order(y.id)]
  y.agg <- df$y
  a.agg <- df$a
  x.agg <- as.matrix(df[,4:ncol(df)])
  offset <- unname(table(y.id))
  
  out2 <- hct_dr(y = y.agg, a = a.agg, x = x.agg, a.vals = a.vals, offset = offset, span.seq = span.seq, k = k)

  est[i,3,] <- out2$estimate
  se[i,2,] <- sqrt(out2$variance)
  
}

cp <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,]) < colMeans(est[,1,]) & (est[i,2,] + 1.96*se[i,]) > colMeans(est[,1,])))

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("TRUE SATE", "DR")
out_cp <- rowMeans(cp, na.rm = T)
names(out_cp) <- a.vals

out_est
out_cp

plot(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "blue", lwd = 2, main = "Exposure Response Curve", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "green", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2)
