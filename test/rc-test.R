### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(parallel)

# Code for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/blp.R")
source("D:/Github/causal-me/hct-dr.R")

# simulation arguments
n.sim <- 100
sig_epe <- sqrt(2)
sig_gps <- 1
sig_berk <- 0
sig_pred <- 0.5
out_scen <- "a"
gps_scen <- "a"
span <- NULL

# gen data arguments
l <- 1000 # c(500, 800)
m <- 100 # c(100, 200)
n <- 2000 # c(1000, 4000)

# dr arguments
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")

# initialize output
est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
se <- array(NA, dim = c(n.sim, 4, length(a.vals)))

for (i in 1:n.sim){
  
  print(i)
  
  # generate data
  dat <- gen_data(l = l, m = m, n = n, 
                  sig_gps = sig_gps, sig_epe = sig_epe,
                  sig_berk = sig_berk, sig_pred = sig_pred, 
                  out_scen = out_scen, gps_scen = gps_scen)
  
  s.id <- dat$s.id
  y.id <- dat$y.id
  y <- dat$y
  s <- dat$s*rbinom(l, 1, 0.1)
  s[s == 0] <- NA
  x <- dat$x
  w <- dat$w
  a <- dat$a
  shat <- dat$shat
  
  # remove clusters w/o exposure data
  ungroup <- which(!(1:m %in% unique(s.id)))
  if (length(ungroup)!= 0) {
    y <- dat$y[!(y.id %in% ungroup)]
    x <- dat$x[!(y.id %in% ungroup),]
    t <- dat$t[!(y.id %in% ungroup)]
    y.id <- dat$y.id[!(y.id %in% ungroup)]
  }
    
  # set up evaluation points & matrices for predictions
  ws <- data.frame(w, shat)
  ws.tmp <- data.frame(ws[!is.na(s),])
  s.tmp <- s[!is.na(s)]
  colnames(ws.tmp) <- colnames(ws) <- c(colnames(w), "expos")
  
  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = s.tmp, X = ws.tmp, SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "SL.gam"))
  stilde <- c(predict(mumod, newdata = ws)$pred)
  
  id <- unique(s.id)[order(unique(s.id))]
  
  # naive w/o pred
  
  zhat <- aggregate(shat, by = list(s.id), mean)[,2]
  zhat_y <- rep(NA, length(y.id))
  
  for (g in id)
    zhat_y[y.id == g] <- zhat[id == g]
  
  naive_wo <- hct_dr(y = y, a = zhat_y, x = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  # naive w/ pred
  
  z <- aggregate(stilde, by = list(s.id), mean)[,2]
  z_y <- rep(NA, length(y.id))
  
  for (g in id)
    z_y[y.id == g] <- zhat[id == g]
  
  naive_w <- hct_dr(y = y, a = z_y, x = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  # blp w/o pred
  
  ahat_y <- blp(s = shat, s.id = s.id, y.id = y.id)$a_y
  blp_wo <- hct_dr(y = y, a = ahat_y, x = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  # blp w/ pred
  
  a_y <- blp(s = stilde, s.id = s.id, y.id = y.id)$a_y
  blp_w <- hct_dr(y = y, a = a_y, x = x, y.id = y.id, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  est[i,1,] <- predict_example(a = a.vals, x = x, y.id = y.id, out_scen = out_scen)
  est[i,2,] <- naive_wo$estimate
  est[i,3,] <- naive_w$estimate
  est[i,4,] <- blp_wo$estimate
  est[i,5,] <- blp_w$estimate

  se[i,1,] <- sqrt(naive_wo$variance)
  se[i,2,] <- sqrt(naive_w$variance)
  se[i,3,] <- sqrt(blp_wo$variance)
  se[i,4,] <- sqrt(blp_w$variance)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("SAMPLE ERC", "NAIVE WO", "NAIVE W", "BLP W", "BLP WO")

cp_naive_wo <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]))

cp_naive_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]))

cp_blp_wo <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]))

cp_blp_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & (est[i,5,] + 1.96*se[i,4,]) > est[i,1,]))

out_cp <- rbind(rowMeans(cp_naive_wo, na.rm = T), rowMeans(cp_naive_wo, na.rm = T),
                rowMeans(cp_blp_wo, na.rm = T), rowMeans(cp_blp_wo, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("NAIVE WO", "NAIVE W", "BLP W", "BLP WO")

out_est
out_cp

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "red", lwd = 2,
      main = "Exposure Response Curve without Berkson Error", xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "purple", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "blue", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[4,], type = "l", col = "yellow", lwd = 2)
lines(a.vals, colMeans(est, na.rm = T)[5,], type = "l", col = "green", lwd = 2)
legend(1.8, 0.2, legend=c("Sample ERC", "Naive W/O Prediction", "Naive W/ Prediction", "BLP W/O Prediction", "BLP W/ Prediction"),
       col=c("red", "purple", "blue", "yellow", "green"), lwd=2, cex=0.8)
