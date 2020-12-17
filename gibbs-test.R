### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)

# Code for generating and fitting data
source("D:/Dropbox (Personal)/ERC-EPE/Code/gen-data.R")
source("D:/Dropbox (Personal)/ERC-EPE/Code/gibbs-dr.R")
source("D:/Dropbox (Personal)/ERC-EPE/Code/hct-dr.R")

# simulation arguments
n.sim <- 100
n.adapt <- 100
n.iter <- 1000
sig_epe <- 1
sig_gps <- 2
prob <- 0.1

# gen data arguments
l <- 2000 # c(500, 800)
m <- 200 # c(100, 200)
n <- 2000 # c(1000, 4000)

# dr arguments
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
trim <- 0.01
thin <- 10

# prior values
a <- 0.001 # gamma shape
b <- 0.001 # gamma rate
c <- 1e5 # normal rate
h <- 1 # mh tuning

est <- array(NA, dim = c(n.sim, 2, length(a.vals)))
se1 <- se2 <- cp <- matrix(NA, n.sim, length(a.vals))

for (i in 1:n.sim){
  
  print(i)
  
  dat <- gen_agg_data(l = l, m = m, n = n,  prob = prob)
  
  s.id <- dat$s.id
  y.id <- dat$y.id
  
  ungroup <- which(!(1:m %in% unique(s.id)))
  
  y <- dat$y
  s <- dat$s
  t <- dat$t
  x <- dat$x
  w <- dat$w
  
  if (length(ungroup)!= 0) {
    y <- dat$y[!(y.id %in% ungroup)]
    x <- dat$x[!(y.id %in% ungroup),]
    y.id <- dat$y.id[!(y.id %in% ungroup)]
  }
    
  fmla.s <- formula("~ w1 + w2 + w3 + w4")
  fmla.z <- formula("~ x1 + x2 + x3 + x4")
  fmla.y <- formula("~ x1 + x2 + x3 + x4 + z*(x1 + x2 + x3 + x4)")
  
  out <- gibbs_dr(s = s, t = t, w = w, s.id = s.id, y = y, x = x, y.id = y.id, 
                  fmla.s = fmla.s, fmla.z = fmla.z, fmla.y = fmla.y, 
                  a = a, b = b, c = c, h = h, a.vals = a.vals, thin = 5)
  
  est[i,1,] <- predict.example(a = a.vals, x = x, id = y.id)
  est[i,2,] <- out$est
  se1[i,] <- out$se1
  se2[i,] <- out$se2
  
  # Check Mixing of MH sampler
  # par(mfrow = c(3, 4))
  # for (i in 1:ncol(out$mcmc$gamma))
  #   plot(out$mcmc$gamma[,i], type = "l",
  #        xlab = "iterations")
  # title("Trace Plots", outer = TRUE)
  
}

cp <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se1[i,]) < est[i,1,] & (est[i,2,] + 1.96*se1[i,]) > est[i,1,]))

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
