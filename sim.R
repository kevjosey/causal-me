### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(parallel)

# Code for generating and fitting data
source("D:/Github/causal-me/gen-data.R")
source("D:/Github/causal-me/gibbs-dr.R")
source("D:/Github/causal-me/simex-dr.R")
source("D:/Github/causal-me/hct-dr.R")
source("D:/Github/causal-me/mclapply-hack.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # gen_data arguments
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  m <- scenario$m
  n <- scenario$n
  prob <- scenario$prob
  sig_agg <- scenario$sig_agg
  sig_gps <- scenario$sig_gps
  sig_pred <- scenario$sig_pred
  
  # mcmc/prior arguments
  shape <- 1e-3 # gamma shape
  rate <- 1e-3 # gamma rate
  scale <- 1e5 # normal scale
  thin <- 20
  n.adapt <- 100
  n.iter <- 1000
  
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
    a <- dat$a
    w <- dat$w
    star <- dat$star
    
    # validation subset
    s <- dat$s*rbinom(m, 1, prob)
    s[s == 0] <- NA
    
    # remove any clusters w/o exposure data (random process)
    keep <- which((id %in% unique(s.id)))
    if (length(keep)!= 0) {
      y <- dat$y[(id %in% keep)]
      x <- dat$x[(id %in% keep),]
      a <- dat$a[(id %in% keep)]
      wts <- dat$wts[(id %in% keep)]
      id <- dat$id[(id %in% keep)]
    }
    
    stilde <- pred(s = s, star = star, w = w, sl.lib = sl.lib)
    
    a_naive_wo <- aggregate(star, by = list(s.id), mean)
    a_naive_w <- aggregate(stilde, by = list(s.id), mean)
    a_wo <- blp(s = stilde, s.id = s.id) 
    a_w <- blp(s = stilde, s.id = s.id)
    
    blp_naive_wo <- erc(y = y, a = a_naive_wo, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
    blp_naive_w <- erc(y = y, a = a_naive_w, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
    blp_wo <- erc(y = y, a = a_wo, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
    blp_w <- erc(y = y, a = a_w, x = x, wts = wts, a.vals = a.vals, sl.lib = sl.lib, span = span)
    
    est[i,1,] <- predict_example(a = a.vals, x = x, id = id, out_scen = out_scen)
    est[i,2,] <- blp_naive_wo$estimate
    est[i,3,] <- blp_naive_w$estimate
    est[i,4,] <- blp_naive_wo$estimate
    est[i,5,] <- colMeans(gibbs_est_x)
    
    se[i,1,] <- sqrt(blp_w$variance)
    se[i,2,] <- sqrt(var_w)
    se[i,3,] <- sqrt(blp_x$variance)
    se[i,4,] <- sqrt(var_x)
    
  }

  out_est <- colMeans(est, na.rm = T)
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("TRUE SATE", "NAIVE WO", "BLP W", "SIMEX", "MI")
  
  out_se <- colMeans(se, na.rm = T)
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("NAIVE", "SI", "SIMEX", "MI")
  
  cp_list <- lapply(1:n.sim, function(i,...)
    list(naive = as.numeric((est[i,2,] - 1.96*se[i,1,]) < est[i,1,] & (est[i,2,] + 1.96*se[i,1,]) > est[i,1,]),
          si = as.numeric((est[i,3,] - 1.96*se[i,2,]) < est[i,1,] & (est[i,3,] + 1.96*se[i,2,]) > est[i,1,]),
          simex = as.numeric((est[i,4,] - 1.96*se[i,3,]) < est[i,1,] & (est[i,4,] + 1.96*se[i,3,]) > est[i,1,]),
          mi = as.numeric((est[i,5,] - 1.96*se[i,4,]) < est[i,1,] & (est[i,5,] + 1.96*se[i,4,]) > est[i,1,])))
  
  cp_naive <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$naive)), ncol = length(cp_list)), na.rm = TRUE)
  cp_si <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$si)), ncol = length(cp_list)), na.rm = TRUE)
  cp_simex <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$simex)), ncol = length(cp_list)), na.rm = TRUE)
  cp_mi <- rowMeans(matrix(unlist(lapply(cp_list, function(x) x$mi)), ncol = length(cp_list)), na.rm = TRUE)
  
  out_cp <- rbind(NAIVE = cp_naive, SI = cp_si, SIMEX = cp_simex, MI = cp_mi)
  colnames(out_cp) <- a.vals

  out <- list(est = out_est, se = out_se, cp = out_cp)
  return(out)
  
}
  
# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(-1, 3, by = 0.25)
sl.lib <- c("SL.mean", "SL.glm", "SL.glm.interaction")
n.sim <- 100

n <- c(2000)
m <- c(100)
l <- c(500)
out_scen <- c("a", "b")
gps_scen <- c("a", "b")

scen_mat <- expand.grid(n = n, m = m, l = l, out_scen = out_scen, gps_scen = gps_scen)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
rslt <- mclapply.hack(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib, mc.cores = 2)
rslt <- list(rslt = rslt, scen_idx = scen_mat)

save(rslt, file = "D:/Github/causal-me/output/rslt.RData")

pdf(file = "D:/Github/causal-me/output/ERC-plot.pdf")
par(mfrow = c(2, 2))
for (k in 1:length(rlist_est)) {
  
  plot(a.vals, rslt[[1]][[k]]$est[1,], type = "l", col = "red", lwd = 2,
       main = paste("scenario", k), xlab = "Exposure", ylab = "Probability of Event", ylim = c(0,1))
  lines(a.vals, rslt[[1]][[k]]$est[2,], type = "l", col = "orange", lwd = 2)
  lines(a.vals, rslt[[1]][[k]]$est[3,], type = "l", col = "green", lwd = 2)
  lines(a.vals, rslt[[1]][[k]]$est[4,], type = "l", col = "blue", lwd = 2)
  lines(a.vals, rslt[[1]][[k]]$est[5,], type = "l", col = "purple", lwd = 2)

  if (k == length(rslt[[1]]))
    legend(1, 0.4, legend=c("Sample ERC", "Naive", "SI", "SIMEX", "MI"), 
           col=c("red", "orange", "green", "blue", "purple"), lwd=2, cex=0.8)
  
}
dev.off()
