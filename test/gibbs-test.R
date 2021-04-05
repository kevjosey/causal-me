### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(gam)
library(parallel)

# Code for generating and fitting data
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/gibbs-sampler.R")
source("~/Github/causal-me/mclapply-hack.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/erc.R")

# simulation arguments
n.sim <- 200
sig_gps <- 1
sig_agg <- sqrt(2)
sig_pred <- sqrt(0.5)
gps_scen <- "a"
out_scen <- "a"
pred_scen <- "b"
span <- 0.5

# gen data arguments
m <- 1000 # c(500, 800)
n <- 200 # c(100, 200)

# gibbs sampler stuff
thin <- 50
n.iter <- 5000
n.adapt <- 1000
h.a <- 1
h.gamma <- 0.3
deg.num <- 2

# dr arguments
a.vals <- seq(6, 10, by = 0.1)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam")
family <- poisson()

# initialize output
est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
se <- array(NA, dim = c(n.sim, 4, length(a.vals)))

for (i in 1:n.sim){
  
  print(i)
  
  # generate data
  dat <- gen_data(m = m, n = n, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
                  pred_scen = pred_scen, out_scen = out_scen, gps_scen = gps_scen)
  
  # zipcode index
  s.id <- dat$s.id
  id <- dat$id
  offset <- log(dat$offset)
  
  # data
  y <- dat$y
  x <- dat$x
  a <- dat$a
  w <- dat$w
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
    offset <- log(dat$offset[(id %in% keep)])
    id <- dat$id[(id %in% keep)]
  }
  
  s_hat <- pred(s = s, star = star, w = w, sl.lib = sl.lib)
  a_w <- blp(s = s_hat, s.id = s.id)
  a_x <- blp(s = s_hat, s.id = s.id, x = x)
  
  # blp w/ pred
  
  blp_w <- erc(y = y, a = a_w, x = x, offset = offset, family = family,
               a.vals = a.vals, sl.lib = sl.lib, span = span)
  blp_x <- erc(y = y, a = a_x, x = x, offset = offset, family = family,
               a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  # Bayesian analysis
  gibbs_w <- gibbs_dr(s = s, star = star,  y = y, offset = offset,
                      s.id = s.id, id = id, w = w, family = family,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin, 
                      h.a = h.a, h.gamma = h.gamma, deg.num = deg.num,
                      a.vals = a.vals, span = span, mc.cores = 4)
  
  gibbs_x <- gibbs_dr(s = s, star = star, y = y, offset = offset,
                      s.id = s.id, id = id, w = w, x = x, family = family,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin, 
                      h.a = 1, h.gamma = 0.3, deg.num = 2,
                      a.vals = a.vals, span = span, mc.cores = 4)
  
  # estimates
  est[i,1,] <- predict_example(a = a.vals, x = x, out_scen = out_scen)
  est[i,2,] <- blp_w$estimate
  est[i,3,] <- gibbs_w$estimate
  est[i,4,] <- blp_x$estimate
  est[i,5,] <- gibbs_x$estimate

  # standard error
  se[i,1,] <- sqrt(blp_w$variance)
  se[i,2,] <- sqrt(gibbs_w$variance)
  se[i,3,] <- sqrt(blp_x$variance)
  se[i,4,] <- sqrt(gibbs_x$variance)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("SAMPLE ERC","BLP W", "Gibbs W", "BLP X", "Gibbs X")

cp_blp_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < colMeans(est[,1,],na.rm = TRUE) & 
               (est[i,2,] + 1.96*se[i,1,]) > colMeans(est[,1,],na.rm = TRUE)))

cp_gibbs_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,3,] - 1.96*se[i,2,]) < colMeans(est[,1,],na.rm = TRUE) & 
               (est[i,3,] + 1.96*se[i,2,]) > colMeans(est[,1,],na.rm = TRUE)))

cp_blp_x <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,4,] - 1.96*se[i,3,]) < colMeans(est[,1,],na.rm = TRUE) & 
               (est[i,4,] + 1.96*se[i,3,]) > colMeans(est[,1,],na.rm = TRUE)))

cp_gibbs_x <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,5,] - 1.96*se[i,4,]) < colMeans(est[,1,],na.rm = TRUE) & 
               (est[i,5,] + 1.96*se[i,4,]) > colMeans(est[,1,],na.rm = TRUE)))

out_cp <- rbind(rowMeans(cp_blp_w, na.rm = T), rowMeans(cp_gibbs_w, na.rm = T),
                rowMeans(cp_blp_x, na.rm = T), rowMeans(cp_gibbs_x, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("BLP W", "Gibbs W", "BLP X", "Gibbs X")

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "darkgreen", lwd = 2,
     main = "Exposure = a, Outcome = a", xlab = "Exposure", ylab = "Rate of Event", 
     ylim = c(0,0.1))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "red", lwd = 2, lty = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "blue", lwd = 2, lty = 2)
lines(a.vals, colMeans(est, na.rm = T)[4,], type = "l", col = "red", lwd = 2, lty = 3)
lines(a.vals, colMeans(est, na.rm = T)[5,], type = "l", col = "blue", lwd = 2, lty = 3)

legend(6, 0.1, legend=c("Sample ERC", "Single Imputation", "Multiple Imputation", "Without Covariates", "With Covariates"),
       col=c("darkgreen", "red", "blue", "black", "black"),
       lty = c(1,1,1,2,3), lwd=2, cex=0.8)
