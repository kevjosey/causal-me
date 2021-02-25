### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(parallel)

# Code for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/gibbs-sampler.R")
source("D:/Github/causal-me/blp.R")
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
span <- 0.8

# mcmc/prior values
shape <- 1e-5 # gamma shape
rate <- 1e-5 # gamma rate
scale <- 1e5 # normal scale
thin <- 10
n.adapt <- 1000
n.iter <- 10000

# initialize output
est <- array(NA, dim = c(n.sim, 4, length(a.vals)))
se <- array(NA, dim = c(n.sim, 3, length(a.vals)))

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
    
  mcmc <- gibbs_dr(s = s, x = x, s.id = s.id, y.id = y.id,
                  shape = shape, rate = rate, scale = scale, 
                  thin = thin, n.iter = n.iter, n.adapt = n.adapt)
  
  bla <- blp(s = s, x = x, s.id = s.id, y.id = y.id)
  
  gibbs_si <- hct_dr(y = y, a = colMeans(mcmc$amat_y), x = x, 
                     y.id = y.id, a.vals = a.vals, span = span, sl.lib = sl.lib)
  
  blp_si <- hct_dr(y = y, a = bla$a_y, x = x, 
                   y.id = y.id, a.vals = a.vals, span = span, sl.lib = sl.lib)
  
  t <- aggregate(s, by = list(s.id), mean)[,2]
  id <- unique(s.id)[order(unique(s.id))]
  t_y <- rep(NA, length(y.id))
  
  for (g in id) {
    t_y[y.id == g] <- t[id == g]
  }
  
  naive <- hct_dr(y = y, a = t_y, x = x, y.id = y.id, a.vals = a.vals, span = span, sl.lib = sl.lib)
  
  est[i,1,] <- predict_example(a = a.vals, x = x, y.id = y.id, out_scen = out_scen)
  est[i,2,] <- naive$estimate
  est[i,3,] <- blp_si$estimate
  est[i,4,] <- gibbs_si$estimate
  
  se[i,1,] <- sqrt(naive$variance)
  se[i,2,] <- sqrt(blp_si$variance)
  se[i,3,] <- sqrt(gibbs_si$variance)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("SAMPLE ERC", "NAIVE", "BLP", "GIBBS")

cp_naive <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]))

cp_blp <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]))

cp_gibbs <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]))

out_cp <- rbind(rowMeans(cp_naive, na.rm = T), rowMeans(cp_blp, na.rm = T), rowMeans(cp_gibbs, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("NAIVE", "BLP", "GIBBS")

out_est
out_cp

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2,
      main = "Exposure Response Curve", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "blue", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "green", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[4,], type = "l", col = "purple", lwd = 2)
legend(2, 0.2, legend=c("Sample ERC", "Naive", "BLP", "Gibbs"), col=c("red", "blue", "green", "purple"), lwd=2, cex=0.8)
