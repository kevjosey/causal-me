### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(earth)
library(splines)
library(parallel)

# Code for generating and fitting data
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/bayes-erc.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/erc.R")

# simulation arguments
n.sim <- 100
sig_gps <- 1
sig_agg <- sqrt(2)
sig_pred <- sqrt(0.5)
gps_scen <- "a"
out_scen <- "a"
pred_scen <- "b"
span <- 0.3

# gen data arguments
mult <- 10
n <- 400 
prob <- 0.2

# gibbs sampler stuff
thin <- 100
n.iter <- 5000
n.adapt <- 500
h.a <- 1
h.gamma <- 0.25
deg.num <- 2

# dr arguments
a.vals <- seq(6, 10, by = 0.04)
sl.lib <- c("SL.mean","SL.glm")
family <- poisson()

# initialize output
est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
se <- array(NA, dim = c(n.sim, 2, length(a.vals)))

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
  s_tilde <- dat$star
  
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
  a_hat <- blp(s = s_hat, s.id = s.id, x = x)
  
  # blp
  blp_hat <- erc(y = y, a = a_hat, x = x, offset = offset, family = family,
               a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num)
  
  # Bayesian analysis
  gibbs_hat <- bayes_erc(s = s, star = s_tilde, y = y, offset = offset,
                         s.id = s.id, id = id, w = w, x = x, family = family,
                         n.iter = n.iter, n.adapt = n.adapt, thin = thin, 
                         h.a = h.a, h.gamma = h.gamma, deg.num = deg.num,
                         a.vals = a.vals, span = span)
  
  # estimates
  est[i,1,] <- predict_example(a = a.vals, x = x, out_scen = out_scen)
  est[i,2,] <- blp_hat$estimate
  est[i,3,] <- gibbs_hat$estimate

  # standard error
  se[i,1,] <- sqrt(blp_hat$variance)
  se[i,2,] <- sqrt(gibbs_hat$variance)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("ERF","SI","Bayes")

cp_blp_x <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & 
               (est[i,2,] + 1.96*se[i,1,]) > est[i,1,],na.rm = TRUE))

cp_gibbs_x <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & 
               (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]))

out_cp <- rbind(rowMeans(cp_blp_x, na.rm = T), rowMeans(cp_gibbs_x, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("SI", "Bayes")

plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "darkgreen", lwd = 2,
     main = "Exposure = a, Outcome = a", xlab = "Exposure", ylab = "Rate of Event", 
     ylim = c(0,0.1))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "red", lwd = 2, lty = 1)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "blue", lwd = 2, lty = 1)

legend(6, 0.1, legend=c("True ERF", "Single Imputation", "Bayesian Approach"),
       col=c("darkgreen", "red", "blue"),
       lty = c(1,1,1), lwd=2, cex=0.8)