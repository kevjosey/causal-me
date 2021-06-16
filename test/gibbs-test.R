### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(splines)
library(scales)
library(parallel)

# Code for generating and fitting data
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/bayes-erc.R")
source("~/Github/causal-me/mi-erc.R")
source("~/Github/causal-me/gps-erc.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/erc.R")
source("~/Github/causal-me/auxiliary.R")

# simulation arguments
n.sim <- 100
sig_gps <- sqrt(4)
sig_agg <- sqrt(2)
sig_pred <- sqrt(0.5)
gps_scen <- "a"
out_scen <- "a"
pred_scen <- "b"
span <- 0.2
mult <- 5
n <- 400
prob <- 0.2

# model arguments
a.vals <- seq(6, 10, by = 0.04)
sl.lib <- c("SL.glm") # for single imputation
family <- poisson()
deg.num <- 2

# mcmc arguments
n.iter <- 10000
n.boot <- 500
n.adapt <- 1000
thin <- 100
h.a <- 1
h.gamma <- c(rep(0.02, 5), rep(0.01, 3))
scale <- 1e6
shape <- rate <- 1e-3

# initialize output
est <- array(NA, dim = c(n.sim, 5, length(a.vals)))
se <- cp <- array(NA, dim = c(n.sim, 4, length(a.vals)))

for (i in 1:n.sim){
  
  print(i)
  
  # generate data
  dat <- gen_data(n = n, mult = mult, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
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
  s_tilde <- star <- dat$star
  
  # validation subset
  s <- dat$s*rbinom(mult*n, 1, prob)
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
  
  s_hat <- pred(s = s, star = s_tilde, w = w, sl.lib = sl.lib)
  a_hat <- blp(s = s_hat, s.id = s.id, x = x)$a
  
  z <- aggregate(s_hat, by = list(s.id), mean)[,2]
  
  # naive
  naive_hat <- erc(y = y, a = z, x = x, offset = offset, family = family,
                   a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num)
  
  # real
  tru_hat <- erc(y = y, a = a, x = x, offset = offset, family = family,
                 a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num)
  
  # blp
  blp_hat <- erc(y = y, a = a_hat, x = x, offset = offset, family = family,
                 a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num)
  
  # Bayesian analysis
  gibbs_hat <- gps_erc(s = s, star = s_tilde, y = y, offset = offset,
                      s.id = s.id, id = id, w = w, x = x, family = family, 
                      deg.num = deg.num, a.vals = a.vals, span = span,
                      scale = scale, shape = shape, rate = rate,
                      h.a = h.a, h.gamma = h.gamma, n.iter = n.iter, 
                      n.adapt = n.adapt, thin = thin)
  
  # estimates
  est[i,1,] <- predict_example(a = a.vals, x = x, out_scen = out_scen)
  est[i,2,] <- tru_hat$estimate
  est[i,3,] <- naive_hat$estimate
  est[i,4,] <- blp_hat$estimate
  est[i,5,] <- gibbs_hat$estimate

  # standard error
  se[i,1,] <- sqrt(tru_hat$variance)
  se[i,2,] <- sqrt(naive_hat$variance)
  se[i,3,] <- sqrt(blp_hat$variance)
  se[i,4,] <- sqrt(gibbs_hat$variance)
  
  # coverage
  cp[i,1,] <- as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & 
                           (est[i,2,] + 1.96*se[i,1,]) > est[i,1,])
  cp[i,2,] <- as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & 
                           (est[i,3,] + 1.96*se[i,2,]) > est[i,1,])
  cp[i,3,] <- as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & 
                           (est[i,4,] + 1.96*se[i,3,]) > est[i,1,])
  cp[i,4,] <- as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & 
                           (est[i,5,] + 1.96*se[i,4,]) > est[i,1,])
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("ERF","DR","Naive","SI","MI")

out_cp <- colMeans(cp, na.rm = T)
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("DR","Naive","SI","MI")

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = hue_pal()(5)[1], lwd = 2,
     main = "Exposure = a, Outcome = a", xlab = "Exposure", ylab = "Rate of Event", 
     ylim = c(0,0.15))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = hue_pal()(5)[2], lwd = 2, lty = 1)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = hue_pal()(5)[3], lwd = 2, lty = 1)
lines(a.vals, colMeans(est, na.rm = T)[4,], type = "l", col = hue_pal()(5)[4], lwd = 2, lty = 1)
lines(a.vals, colMeans(est, na.rm = T)[5,], type = "l", col = hue_pal()(5)[5], lwd = 2, lty = 1)

legend(6, 0.15, legend=c("True ERF", "Observed", "Naive", "Single Imputation", "Multiple Imputation"),
       col=hue_pal()(5),
       lty = c(1,1,1,1,1), lwd=2, cex=0.8)
