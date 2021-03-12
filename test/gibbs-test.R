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

# simulation arguments
n.sim <- 100
sig_agg <- sqrt(2)
sig_gps <- 1
sig_berk <- 0
sig_pred <- 1
out_scen <- "a"
gps_scen <- "a"
span <- 0.8

# gen data arguments
l <- 1000 # c(500, 800)
m <- 100 # c(100, 200)
n <- 2000 # c(1000, 4000)

# gibbs sampler stuff
thin = 100
n.iter = 1000
n.adapt = 100

# dr arguments
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")

# initialize output
est <- array(NA, dim = c(n.sim, 4, length(a.vals)))
se <- array(NA, dim = c(n.sim, 3, length(a.vals)))

for (i in 1:n.sim){
  
  print(i)
  
  # generate data
  dat <- gen_data(l = l, m = m, n = n, 
                  sig_gps = sig_gps, sig_agg = sig_agg,
                  sig_berk = sig_berk, sig_pred = sig_pred, 
                  out_scen = out_scen, gps_scen = gps_scen)
  
  # zipcode index
  s.id <- dat$s.id
  y.id <- dat$y.id
  id <- unique(s.id)[order(unique(s.id))]
  
  # data
  y <- dat$y
  x <- dat$x
  w <- dat$w
  a <- dat$a
  shat <- dat$shat
  
  # validation subset
  s <- dat$s*rbinom(l, 1, 0.1)
  s[s == 0] <- NA
  
  # remove clusters w/o exposure data
  ungroup <- which(!(1:m %in% unique(s.id)))
  if (length(ungroup)!= 0) {
    y <- dat$y[!(y.id %in% ungroup)]
    x <- dat$x[!(y.id %in% ungroup),]
    t <- dat$t[!(y.id %in% ungroup)]
    y.id <- dat$y.id[!(y.id %in% ungroup)]
  }
  
  stilde <- pred(s = s, shat = shat, w = w)
  
  ahat <- blp(s = shat, s.id = s.id)$a
  a <- blp(s = stilde, s.id = s.id)$a
  ahat_y <- a_y <- rep(NA, length(y.id))
  
  for (g in id){
    
    ahat_y[y.id == g] <- ahat[id == g]
    a_y[y.id == g] <- a[id == g]
    
  }
  
  gibbs_w <- gibbs_dr(s = shat, s.id = s.id, y.id = y.id, x = x, w = w, t = s,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  out <- apply(gibbs_w$amat_y, 1, hct_dr, y = y, x = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, span = span)
  gibbs_est <- do.call(rbind, lapply(out, function(o) o$estimate))
  
  # blp w/o pred
  
  blp_wo <- hct_dr(y = y, a = ahat_y, x = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  # blp w/ pred
  
  blp_w <- hct_dr(y = y, a = a_y, x = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  est[i,1,] <- predict_example(a = a.vals, x = x, y.id = y.id, out_scen = out_scen)
  est[i,2,] <- colMeans(gibbs_est)
  est[i,3,] <- blp_wo$estimate
  est[i,4,] <- blp_w$estimate
  
  se[i,1,] <- sqrt(apply(gibbs_est, 2, var))
  se[i,2,] <- sqrt(blp_wo$variance)
  se[i,3,] <- sqrt(blp_w$variance)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("SAMPLE ERC", "Gibbs", "BLP W", "BLP WO")

cp_gibbs <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]))

cp_blp_wo <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]))

cp_blp_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & (est[i,5,] + 1.96*se[i,4,]) > est[i,1,]))

out_cp <- rbind(rowMeans(cp_gibbs, na.rm = T), rowMeans(cp_blp_wo, na.rm = T), rowMeans(cp_blp_w, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("NAIVE WO", "NAIVE W", "BLP W", "BLP WO")

out_est
out_cp

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2,
     main = "Exposure Response Curve without Berkson Error", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "purple", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "blue", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[4,], type = "l", col = "darkgreen", lwd = 2)
legend(1.8, 0.2, legend=c("Sample ERC", "Gibbs", "BLP W/O Prediction", "BLP W/ Prediction"),
       col=c("red", "purple", "blue", "darkgreen"), lwd=2, cex=0.8)
