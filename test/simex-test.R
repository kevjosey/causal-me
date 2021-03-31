### Test DR estimator

rm(list = ls())

## Preliminaries

library(mvtnorm)
library(cobalt)
library(SuperLearner)
library(data.table)
library(parallel)

# Code for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/simex.R")
source("D:/Github/causal-me/blp.R")
source("D:/Github/causal-me/gibbs-sampler.R")
source("D:/Github/causal-me/erc.R")
source("D:/Github/causal-me/mclapply-hack.R")

n.sim <- 100
n.boot <- 50

# dimensions
m <- 500 # c(100, 200)
n <- 100 # c(1000, 4000)

# simulation arguments
n.sim <- 200
sig_gps <- 1
sig_agg <- sqrt(2)
sig_pred <- sqrt(0.5)
gps_scen <- "a"
out_scen <- "a"
pred_scen <- "b"

# other preliminaries
span <- 0.5
degree <- 2
lambda <- seq(0.1, 2.1, by = 0.25)

# gibbs sampler stuff
thin <- 20
n.iter <- 1000
n.adapt <- 100

a.vals <- seq(6, 10, by = 0.5)
sl.lib <- c("SL.mean","SL.glm","SL.glm.interaction")
family <- poisson()

est <- array(NA, dim = c(n.sim, 5, length(a.vals)))

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
  
  z_hat_tmp <- aggregate(s_hat, by = list(s.id), mean)
  z_hat <- rep(NA, length(id))
  
  for (g in id) {
    
    z_hat[id == g] <- z_hat_tmp[z_hat_tmp[,1] == g,2]
    
  }
  
  naive <- erc(y = y, a = z_hat, x = x, offset = offset, family = family,
               a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  
  a_hat <- blp(s = s_hat, s.id = s.id, x = x)
  
  blp_hat <- erc(y = y, a = a_hat, x = x, offset = offset, family = family,
                 a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  w_new <- model.matrix(~ 0 + star*w)

  gibbs_hat <- gibbs_dr(s = s, star = star, s.id = s.id, id = id, w = w_new, x = x,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  
  out_gibbs <- mclapply.hack(1:nrow(gibbs_hat$amat), function(k, ...){
    
    erc(y = y, a = gibbs_hat$amat[k,], x = x, offset = offset, family = family,
        a.vals = a.vals, sl.lib = sl.lib, span = span)
    
  }, mc.cores = 4)
  
  out <- simex(s = s_hat, y = y, x = x, id = id, s.id = s.id, family = family,
               offset = offset, a.vals = a.vals, n.boot = n.boot,
               degree = degree, lambda = lambda, span = span, mc.cores = 4)
  
  gibbs_est <- do.call(rbind, lapply(out_gibbs, function(o) o$estimate))

  est[i,1,] <- predict_example(a.vals = a.vals, x = x, out_scen = "a")
  est[i,2,] <- naive$estimate
  est[i,3,] <- blp_hat$estimate
  est[i,4,] <- colMeans(gibbs_est)
  est[i,5,] <- out$estimate
  
}

plot(a.vals, colMeans(est[,1,], na.rm = T), type = "l", col = "darkgreen", lwd = 2, main = "Exposure Response Curve",
     xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,0.1))
lines(a.vals, colMeans(est[,2,], na.rm = T), type = "l", col = "red", lwd = 2)
lines(a.vals, colMeans(est[,3,], na.rm = T), type = "l", col = "blue", lwd = 2)
lines(a.vals, colMeans(est[,4,], na.rm = T), type = "l", col = "orange", lwd = 2)
lines(a.vals, colMeans(est[,5,], na.rm = T), type = "l", col = "purple", lwd = 2)

legend(6, 0.1, legend=c("Sample ERC", "Naive", "BLP", "MI", "SIMEX"),
       col=c("darkgreen", "red", "blue", "black", "black"),
       lty = c(1,1,1,2,3), lwd=2, cex=0.8)
