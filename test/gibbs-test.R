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
source("D:/Github/causal-me/mclapply-hack.R")
source("D:/Github/causal-me/blp.R")
source("D:/Github/causal-me/erc.R")

# simulation arguments
n.sim <- 100
sig_gps <- 1
sig_agg <- sqrt(2)
sig_pred <- sqrt(0.5)
out_scen <- "a"
gps_scen <- "a"
span <- NULL

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
est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
se <- array(NA, dim = c(n.sim, 2, length(a.vals)))

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
    wts <- dat$wts[(id %in% keep)]
    id <- dat$id[(id %in% keep)]
  }
  
  stilde <- pred(s = s, star = star, w = w, sl.lib = sl.lib)
  a_w <- blp(s = stilde, s.id = s.id)
  
  w_new <- model.matrix(~ 0 + star*w)
  gibbs_w <- gibbs_dr(s = s, star = star, s.id = s.id, id = id, w = w_new,
                      n.iter = n.iter, n.adapt = n.adapt, thin = thin)

  # blp w/ pred
  
  blp_w <- erc(y = y, a = a_w, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)

  # gibbs w/ pred
  
  a_list <- split(gibbs_w$amat, seq(nrow(gibbs_w$amat)))
  
  out_w <- mclapply.hack(1:length(a_list), function(k, ...){
    
    erc(y = y, a = a_list[[k]], x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
    
  }, mc.cores = 3)
  
  gibbs_est_w <- do.call(rbind, lapply(out_w, function(o) o$estimate))
  gibbs_var_w <- do.call(rbind, lapply(out_w, function(o) o$variance))
  
  # estimates
  
  est[i,1,] <- predict_example(a = a.vals, x = x, id = id, out_scen = out_scen)
  est[i,2,] <- blp_w$estimate
  est[i,3,] <- colMeans(gibbs_est_w)

  #standard error
  
  var_w <- colMeans(gibbs_var_w) + (1 + 1/nrow(gibbs_est_w))*apply(gibbs_est_w, 2, var)

  se[i,1,] <- sqrt(blp_w$variance)
  se[i,2,] <- sqrt(var_w)
  
}

out_est <- colMeans(est, na.rm = T)
colnames(out_est) <- a.vals
rownames(out_est) <- c("SAMPLE ERC", "BLP", "Gibbs")

cp_blp_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]))

cp_gibbs_w <- sapply(1:n.sim, function(i,...)
  as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]))

out_cp <- rbind(rowMeans(cp_blp_w, na.rm = T), rowMeans(cp_gibbs_w, na.rm = T))
colnames(out_cp) <- a.vals
rownames(out_cp) <- c("BLP", "Gibbs")

save(out_est, file = "D:/Github/causal-me/output/dr_rslt_a_a.RData")
save(out_cp, file = "D:/Github/causal-me/output/cp_rslt_a_a.RData")

pdf("D:/Github/causal-me/output/rslt_a_a.pdf")
plot(a.vals, colMeans(est, na.rm = T)[1,], type = "l", col = "darkgreen", lwd = 2,
     main = "Exposure = a, Outcome = a", xlab = "Exposure", ylab = "Rate of Event", ylim = c(0.1,0.5))
lines(a.vals, colMeans(est, na.rm = T)[2,], type = "l", col = "red", lwd = 2, lty = 2)
lines(a.vals, colMeans(est, na.rm = T)[3,], type = "l", col = "blue", lwd = 2, lty = 2)

legend(-1, 0.5, legend=c("Sample ERC", "Single Imputation", "Multiple Imputation"),
       col=c("darkgreen", "red", "blue"),
       lty = c(1,2,2), lwd=2, cex=0.8)
dev.off()
