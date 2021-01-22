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
degree <- 3
lambda <- seq(0.1, 2.1, by = 0.25)

# initialize output
est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
se <- cp <- matrix(NA, n.sim, length(a.vals))

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
    a_y <- dat$a_y[!(y.id %in% ungroup)]
    y.id <- dat$y.id[!(y.id %in% ungroup)]
  }
    
  fmla.s <- formula("~ w1 + w2 + w3 + w4")
  fmla.a <- formula("~ x1 + x2 + x3 + x4")
  
  out <- gibbs_dr(s = s, y = y, x = x, s.id = s.id, y.id = y.id, 
                  fmla.s = fmla.s, fmla.a = fmla.a, sl.lib = sl.lib,
                  shape = shape, rate = rate, scale = scale, thin = thin,
                  n.iter = n.iter, n.adapt = n.adapt, n.boot = n.boot, 
                  span = span, lambda = lambda, a.vals = a.vals, degree = degree)
  
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
