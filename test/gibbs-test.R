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
n.sim <- 1000
sig_gps <- 1
sig_agg <- sqrt(2)
sig_pred <- sqrt(0.5)
gps_scen <- "a"
out_scen <- "a"
span <- 0.5

# gen data arguments
m <- 1000 # c(500, 800)
n <- 100 # c(100, 200)

# gibbs sampler stuff
thin = 20
n.iter = 1000
n.adapt = 100

# dr arguments
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam")
formula <- as.formula("y ~ s(x1, x2, x3, x4, df = 3)")

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
  offset <- log(dat$wts)
  
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
    offset <- log(dat$wts[(id %in% keep)])
    id <- dat$id[(id %in% keep)]
  }
  
  stilde <- pred(s = s, star = star, w = w, sl.lib = sl.lib)
  a_w <- blp(s = stilde, s.id = s.id)
  a_x <- blp(s = stilde, s.id = s.id, x = x)
  
  w_new <- model.matrix(~ 0 + star*w)
  gibbs_w <- gibbs_dr(s = s, star = star, s.id = s.id, id = id, w = w_new,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  
  gibbs_x <- gibbs_dr(s = s, star = star, s.id = s.id, id = id, w = w_new, x = x,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin)

  # blp w/ pred
  
  blp_w <- erc(y = y, a = a_w, x = x, offset = offset, a.vals = a.vals, sl.lib = sl.lib, span = span)
  blp_x <- erc(y = y, a = a_x, x = x, offset = offset, a.vals = a.vals, sl.lib = sl.lib, span = span)
  
  # gibbs w/ pred
  
  a_list_w <- split(gibbs_w$amat, seq(nrow(gibbs_w$amat)))
  a_list_x <- split(gibbs_x$amat, seq(nrow(gibbs_x$amat)))
  
  out_w <- mclapply.hack(1:length(a_list_x), function(k, ...){
    
    erc(y = y, a = a_list_w[[k]], x = x, offset = offset, a.vals = a.vals, sl.lib = sl.lib, span = span)
    
  }, mc.cores = 3)
  
  out_x <- mclapply.hack(1:length(a_list_x), function(k, ...){
    
    erc(y = y, a = a_list_x[[k]], x = x, offset = offset, a.vals = a.vals, sl.lib = sl.lib, span = span)
    
  }, mc.cores = 3)
  
  gibbs_est_w <- do.call(rbind, lapply(out_w, function(o) o$estimate))
  gibbs_var_w <- do.call(rbind, lapply(out_w, function(o) o$variance))
  
  gibbs_est_x <- do.call(rbind, lapply(out_x, function(o) o$estimate))
  gibbs_var_x <- do.call(rbind, lapply(out_x, function(o) o$variance))
  
  # estimates
  
  est[i,1,] <- predict_example(a = a.vals, x = x, id = id, out_scen = out_scen)
  est[i,2,] <- blp_w$estimate
  est[i,3,] <- colMeans(gibbs_est_w)
  est[i,4,] <- blp_x$estimate
  est[i,5,] <- colMeans(gibbs_est_x)

  #standard error
  
  var_w <- colMeans(gibbs_var_w) + (1 + 1/nrow(gibbs_est_w))*apply(gibbs_est_w, 2, var)
  var_x <- colMeans(gibbs_var_x) + (1 + 1/nrow(gibbs_est_x))*apply(gibbs_est_x, 2, var)
  
  se[i,1,] <- sqrt(blp_w$variance)
  se[i,2,] <- sqrt(var_w)
  se[i,3,] <- sqrt(blp_x$variance)
  se[i,4,] <- sqrt(var_x)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("SAMPLE ERC","BLP W", "Gibbs W", "BLP X", "Gibbs X")

cp_blp_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]))

cp_gibbs_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]))

cp_blp_x <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]))

cp_gibbs_x <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & (est[i,5,] + 1.96*se[i,4,]) > est[i,1,]))

out_cp <- rbind(rowMeans(cp_blp_w, na.rm = T), rowMeans(cp_gibbs_w, na.rm = T),
                rowMeans(cp_blp_x, na.rm = T), rowMeans(cp_gibbs_x, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("BLP W", "Gibbs W", "BLP X", "Gibbs X")

save(out_est, file = "~/Github/causal-me/output/dr_rslt_a_a.RData")
save(out_cp, file = "~/Github/causal-me/output/cp_rslt_a_a.RData")

pdf("~/Github/causal-me/output/rslt_a_a.pdf")
plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "darkgreen", lwd = 2,
     main = "Exposure = a, Outcome = a", xlab = "Exposure", ylab = "Rate of Event", ylim = c(0.1,0.5))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "red", lwd = 2, lty = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "blue", lwd = 2, lty = 2)
lines(a.vals, colMeans(est, na.rm = T)[4,], type = "l", col = "red", lwd = 2, lty = 3)
lines(a.vals, colMeans(est, na.rm = T)[5,], type = "l", col = "blue", lwd = 2, lty = 3)

legend(-1, 0.5, legend=c("Sample ERC", "Single Imputation", "Multiple Imputation", "Without Covariates", "With Covariates"),
       col=c("darkgreen", "red", "blue", "black", "black"),
       lty = c(1,1,1,2,3), lwd=2, cex=0.8)
dev.off()
