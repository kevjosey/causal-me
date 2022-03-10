### Test DR estimator

rm(list = ls())

## Preliminaries
library(mvtnorm)
library(SuperLearner)
library(scales)
library(parallel)
library(dbarts)
library(abind)

# Code for generating and fitting data
source("~/Github/causal-me/sim/gen-data.R")
source("~/Github/causal-me/erf.R")
source("~/Github/causal-me/bart-erf.R")
source("~/Github/causal-me/bayes-erf.R")
source("~/Github/causal-me/auxiliary.R")

# simulation arguments
n.sim <- 100
sig_gps <- 2
sig_agg <- sqrt(2)
sig_pred <- 1
gps_scen <- "a"
out_scen <- "b"
pred_scen <- "a"
span <- 0.1
mult <- 5
n <- 400
prob <- 0.1

# model arguments
a.vals <- seq(6, 14, by = 0.04)
family <- poisson()

# mcmc arguments
n.iter <- 2000
n.adapt <- 2000
thin <- 20
h.a <- 0.5
scale <- 1e6
shape <- rate <- 1e-3

start <- Sys.time()

out <- mclapply(1:n.sim, function(i, ...){
  
  dat <- gen_data(n = n, mult = mult, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
                  pred_scen = pred_scen, out_scen = out_scen, gps_scen = gps_scen)

  # zipcode index
  s.id <- dat$s.id
  id <- dat$id
  weights <- dat$offset
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
  
  # exposure predictions
  s_hat <- pred(s = s, star = s_tilde, w = w, sl.lib = "SL.glm")
  z_hat <- aggregate(s_hat, by = list(s.id), mean)[,2]
  z_tilde <- aggregate(s_tilde, by = list(s.id), mean)[,2]
  
  # naive
  naive_hat <- try(erf(y = y, a = z_tilde, x = x, offset = offset, weights = weights, 
                       family = family, a.vals = a.vals, span = span,
                       n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
  
  # real
  rc_hat <- try(erf(y = y, a = z_hat, x = x, offset = offset, weights = weights,
                    family = family, a.vals = a.vals, span = span,
                    n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
  
  # BART Approach
  bart_hat <- try(bart_erf(s = s, t = s_tilde, y = y, offset = offset, weights = weights,
                           s.id = s.id, id = id, w = w, x = x, family = family,
                           a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                           h.a = h.a, n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
  
  # Bayes DR Approach
  bayes_hat <- try(bayes_erf(s = s, t = s_tilde, y = y, offset = offset, weights = weights,
                             s.id = s.id, id = id, w = w, x = x, family = family, df = 4,
                             a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                             h.a = h.a, n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
  
  # estimates
  est <- rbind(predict_example(a = a.vals, x = x, out_scen = out_scen),
               if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))},
               if (!inherits(rc_hat, "try-error")) {rc_hat$estimate} else {rep(NA, length(a.vals))},
               if (!inherits(bart_hat, "try-error")) {bart_hat$tree_estimate} else {rep(NA, length(a.vals))},
               if (!inherits(bart_hat, "try-error")) {bart_hat$smooth_estimate} else {rep(NA, length(a.vals))},
               if (!inherits(bayes_hat, "try-error")) {bayes_hat$dr_estimate} else {rep(NA, length(a.vals))})
  
  #standard error
  se <- rbind(if (!inherits(naive_hat, "try-error")) {sqrt(naive_hat$variance)} else {rep(NA, length(a.vals))},
              if (!inherits(rc_hat, "try-error")) {sqrt(rc_hat$variance)} else {rep(NA, length(a.vals))},
              if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$tree_variance)} else {rep(NA, length(a.vals))},
              if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$smooth_variance)} else {rep(NA, length(a.vals))},
              if (!inherits(bayes_hat, "try-error")) {sqrt(bayes_hat$dr_variance)} else {rep(NA, length(a.vals))})
  
  return(list(est = est, se = se))
  
}, mc.cores = 25, mc.preschedule = TRUE)

stop <- Sys.time()
stop - start

est <- abind(lapply(out, function(lst, ...) lst$est), along = 3)
se <- abind(lapply(out, function(lst, ...) lst$se), along = 3)
mu.mat <- matrix(rep(rowMeans(est[1,,]), n.sim), nrow = length(a.vals), ncol = n.sim)

# coverage probability
cp <- list(as.matrix((est[2,,] - 1.96*se[1,,]) < mu.mat & (est[2,,] + 1.96*se[1,,]) > mu.mat),
           as.matrix((est[3,,] - 1.96*se[2,,]) < mu.mat & (est[3,,] + 1.96*se[2,,]) > mu.mat),
           as.matrix((est[4,,] - 1.96*se[3,,]) < mu.mat & (est[5,,] + 1.96*se[3,,]) > mu.mat),
           as.matrix((est[5,,] - 1.96*se[4,,]) < mu.mat & (est[5,,] + 1.96*se[4,,]) > mu.mat),
           as.matrix((est[6,,] - 1.96*se[5,,]) < mu.mat & (est[6,,] + 1.96*se[5,,]) > mu.mat))

out_est <- t(apply(est, 1, rowMeans, na.rm = T))
colnames(out_est) <- a.vals
rownames(out_est) <- c("ERF","NAIVE","RC","BART","LOESS","DR")

out_bias <- t(apply(est[2:6,,], 1, function(x) rowMeans(abs(x - mu.mat), na.rm = T)))
colnames(out_bias) <- a.vals
rownames(out_bias) <- c("NAIVE","RC","BART","LOESS","DR")

out_mse <- t(apply(est[2:6,,], 1, function(x) rowMeans((x - mu.mat)^2, na.rm = T)))
colnames(out_mse) <- a.vals
rownames(out_mse) <- c("NAIVE","RC","BART","LOESS","DR")

out_cp <- do.call(rbind, lapply(cp, rowMeans, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("NAIVE","RC","BART","LOESS", "DR")

plot(a.vals, out_est[1,], type = "l", col = hue_pal()(6)[1], lwd = 2,
     main = "Exposure = b, Outcome = a", xlab = "Exposure", ylab = "Rate of Event", 
     ylim = c(0,0.1))
lines(a.vals, out_est[2,], type = "l", col = hue_pal()(6)[2], lwd = 2, lty = 1)
lines(a.vals, out_est[3,], type = "l", col = hue_pal()(6)[3], lwd = 2, lty = 1)
lines(a.vals, out_est[4,], type = "l", col = hue_pal()(6)[4], lwd = 2, lty = 1)
lines(a.vals, out_est[5,], type = "l", col = hue_pal()(6)[5], lwd = 2, lty = 1)
lines(a.vals, out_est[6,], type = "l", col = hue_pal()(6)[6], lwd = 2, lty = 1)
legend(6, 0.1, legend=c("True ERF","NAIVE","RC","BART","LOESS","DR"),
       col=hue_pal()(6),lty = c(1,1,1,1,1,1), lwd=2, cex=0.8)
