### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(parallel)

# Code for generating and fitting data
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/gibbs-sampler.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/hct-dr.R")
source("~/Github/causal-me/erc.R")

# simulation arguments
n.sim <- 100
sig_agg <- sqrt(2)
sig_gps <- 1
sig_pred <- sqrt(0.5)
out_scen <- "a"
gps_scen <- "a"
span <- 0.75

# gen data arguments
m <- 1000 # c(500, 800)
n <- 100 # c(100, 200)

# gibbs sampler stuff
thin = 20
n.iter = 1000
n.adapt = 100

# dr arguments
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")

# initialize output
est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
se <- array(NA, dim = c(n.sim, 4, length(a.vals)))

for (i in 1:n.sim){
  
  print(i)
  
  # generate data
  dat <- gen_data(m = m, n = n, sig_gps = sig_gps, sig_agg = sig_agg,
                  sig_pred = sig_pred, out_scen = out_scen, gps_scen = gps_scen)
  
  # zipcode index
  s.id <- dat$s.id
  id <- dat$id
  wts <- dat$wts
  
  # data
  y <- dat$y
  x <- dat$x
  w <- dat$w
  a <- dat$a
  star <- dat$star
  
  # validation subset
  s <- dat$s*rbinom(m, 1, 0.1)
  s[s == 0] <- NA
  
  # remove clusters w/o exposure data
  keep <- which((id %in% unique(s.id)))
  if (length(keep)!= 0) {
    y <- dat$y[(id %in% keep)]
    x <- dat$x[(id %in% keep),]
    a <- dat$a[(id %in% keep)]
    wts <- dat$wts[(id %in% keep)]
    id <- dat$id[(id %in% keep)]
  }
  
  stilde <- pred(s = s, star = star, w = w)
  
  ahat <- blp(s = star, s.id = s.id)
  a2 <- blp(s = stilde, s.id = s.id)
  # gibbs_wo <- gibbs_dr(s = star, s.id = s.id, id = id, x = x,
  #                      n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  gibbs_w <- gibbs_dr(s = star, s.id = s.id, id = id, x = x, w = w, t = s,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  
  # blp w/o pred
  
  blp_wo <- erc(y = y, a = ahat, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  # blp w/ pred
  
  blp_w <- erc(y = y, a = a2, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)

  # gibbs w/o pred
  
  # out_wo <- apply(gibbs_wo$amat, 1, erc, y = y, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
  # gibbs_est_wo <- do.call(rbind, lapply(out_wo, function(o) o$estimate))
  
  # gibbs w/ pred
  
  out_w <- apply(gibbs_w$amat, 1, erc, y = y, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
  gibbs_est_w <- do.call(rbind, lapply(out_w, function(o) o$estimate))
  
  est[i,1,] <- predict_example(a = a.vals, x = x, id = id, out_scen = out_scen)
  est[i,2,] <- NA
  est[i,3,] <- colMeans(gibbs_est_w)
  est[i,4,] <- blp_wo$estimate
  est[i,5,] <- blp_w$estimate
  
  se[i,1,] <- NA
  se[i,2,] <- sqrt(apply(gibbs_est_w, 2, var))
  se[i,3,] <- sqrt(blp_wo$variance)
  se[i,4,] <- sqrt(blp_w$variance)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("SAMPLE ERC", "GIBBS WO", "GIBBS W", "BLP WO", "BLP W")

cp_gibbs_wo <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]))

cp_gibbs_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]))

cp_blp_wo <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]))

cp_blp_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & (est[i,5,] + 1.96*se[i,4,]) > est[i,1,]))

out_cp <- rbind(rowMeans(cp_gibbs, na.rm = T), rowMeans(cp_blp_wo, na.rm = T), rowMeans(cp_blp_w, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("GIBBS WO", "GIBBS W", "BLP W", "BLP WO")

out_est
out_cp

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2,
     main = "Exposure Response Curve without Berkson Error", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,0.5))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "purple", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "blue", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[4,], type = "l", col = "darkgreen", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[5,], type = "l", col = "green", lwd = 2)
legend(1.6, 0.1, legend=c("Sample ERC", "Gibbs W/O Prediction", "Gibbs W/ Prediction", 
                          "BLP W/O Prediction", "BLP W/ Prediction"),
       col=c("red", "purple", "blue", "darkgreen", "green"), lwd=2, cex=0.8)
