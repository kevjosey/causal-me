### Test DR estimator

rm(list = ls())

## Preliminaries

library(data.table)
library(mvtnorm)
library(SuperLearner)
library(earth)
library(splines)
library(parallel)
library(abind)

# Code for generating and fitting data
source("~/Github/causal-me/mclapply-hack.R")
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/bayes-erc.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/erc.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # simulation arguments
  sig_gps <- 1
  sig_agg <- sqrt(2)
  sig_pred <- sqrt(0.5)
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  pred_scen <- "b"
  prob <- 0.2
  
  # gen data arguments
  n <- scenario$n # c(500, 800)
  mult <- scenario$mult # c(100, 200)
  
  # gibbs sampler stuff
  thin <- 100
  n.iter <- 5000
  n.adapt <- 500
  h.a <- 1
  h.gamma <- 0.25
  deg.num <- 2
  
  # dr arguments
  span <- ifelse(n == 800, 0.125, 0.25)
  family <- poisson()
  
  # initialize output
  est <- array(NA, dim = c(n.sim, 3, length(a.vals)))
  se <- array(NA, dim = c(n.sim, 2, length(a.vals)))
  
  print(scenario)
  
  # generate data
  dat_list <- replicate(n.sim, gen_data(n = n, mult = mult, sig_gps = sig_gps, sig_agg = sig_agg, sig_pred = sig_pred,
                                        pred_scen = pred_scen, out_scen = out_scen, gps_scen = gps_scen))
  
  out <- mclapply.hack(1:n.sim, function(i,...){

    dat <- dat_list[,i]
    
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
    
    # exposure predictions
    s_hat <- pred(s = s, star = s_tilde, w = w, sl.lib = sl.lib)
    a_hat <- blp(s = s_hat, s.id = s.id, x = x)
    
    # blp w/ pred
    blp_hat <- try(erc(y = y, a = a_hat, x = x, offset = offset, family = family,
                     a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    
    # Bayesian Approach
    gibbs_hat <- try(bayes_erc(s = s, star = s_tilde, y = y, offset = offset,
                               s.id = s.id, id = id, w = w, x = x, family = family,
                               n.iter = n.iter, n.adapt = n.adapt, thin = thin, 
                               h.a = h.a, h.gamma = h.gamma, deg.num = deg.num,
                               a.vals = a.vals, span = span), silent = TRUE)
    
    # estimates
    est <- rbind(predict_example(a = a.vals, x = x, out_scen = out_scen),
                 if (!inherits(blp_hat, "try-error")) {blp_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(gibbs_hat, "try-error")) {gibbs_hat$estimate} else {rep(NA, length(a.vals))})
    
    #standard error
    se <- rbind(if (!inherits(blp_hat, "try-error")) {sqrt(blp_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(gibbs_hat, "try-error")) {sqrt(gibbs_hat$variance)} else {rep(NA, length(a.vals))})
   
    return(list(est = est, se = se))
     
  }, mc.cores = 8)
  
  est <- abind(lapply(out, function(lst, ...) lst$est), along = 3)
  se <- abind(lapply(out, function(lst, ...) lst$se), along = 3)
  
  out_est <- t(apply(est, 1, rowMeans, na.rm = T))
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERF", "SI", "Bayes")
  
  compare <- matrix(rowMeans(est[1,,]), nrow = length(a.vals), ncol = n.sim)
  out_bias <- t(apply(est[2:3,,], 1, function(x) rowMeans(abs(x - compare), na.rm = T)))
  colnames(out_bias) <- a.vals
  rownames(out_bias) <- c("SI", "Bayes")
  
  out_sd <- t(apply(est[2:3,,], 1, function(x) apply(x, 1, sd, na.rm = T)))
  colnames(out_sd) <- a.vals
  rownames(out_sd) <- c("SI", "Bayes")
  
  out_se <- t(apply(se, 1, rowMeans, na.rm = T))
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("SI", "Bayes")
  
  cp_blp_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[2,,i] - 1.96*se[1,,i]) < rowMeans(compare) & 
                 (est[2,,i] + 1.96*se[1,,i]) > rowMeans(compare)))
  
  cp_gibbs_x <- sapply(1:n.sim, function(i,...)
    as.numeric((est[3,,i] - 1.96*se[2,,i]) < rowMeans(compare) & 
                 (est[3,,i] + 1.96*se[2,,i]) > rowMeans(compare)))
  
  out_cp <- rbind(rowMeans(cp_blp_x, na.rm = T), rowMeans(cp_gibbs_x, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_cp) <- c("SI", "Bayes")
  
  return(list(est = out_est, bias = out_bias, sd = out_sd, se = out_se, cp = out_cp))
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.04)
sl.lib <- c("SL.mean","SL.glm")
n.sim <- 1000

n <- c(400,800)
mult <- c(5,10)
gps_scen <- c("a", "b")
out_scen <- c("a", "b")

scen_mat <- expand.grid(n = n, mult = mult, gps_scen = gps_scen, out_scen = out_scen)
scenarios <- lapply(seq_len(nrow(scen_mat)), function(i) scen_mat[i,])
est <- lapply(scenarios, simulate, n.sim = n.sim, a.vals = a.vals, sl.lib = sl.lib)
rslt <- list(est = est, scen_idx = scen_mat)

save(rslt, file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/sim_2.RData")

# Summary Plot

load(file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/sim_2.RData")
plotnames <- c("GPS: \"a\"; Outcome: \"a\"",
               "GPS: \"b\"; Outcome: \"a\"",
               "GPS: \"a\"; Outcome: \"b\"",
               "GPS: \"b\"; Outcome: \"b\"")
idx <- c(4,8,12,16)

filename <- paste0("~/Dropbox (Personal)/Projects/ERC-EPE/Output/plot_2.pdf")
pdf(file = filename, width = 10, height = 10)
par(mfrow = c(2,2))

for (k in 1:4){
  
  plot(a.vals, rslt$est[[idx[k]]]$est[1,], type = "l", col = "darkgreen", lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[k],
       ylim = c(0,0.08))
  grid(lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[2,], type = "l", col = "red", lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[3,], type = "l", col = "blue", lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[2,] -
          1.96*rslt$est[[idx[k]]]$se[1,], type = "l", col = "red", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[idx[k]]]$est[2,] + 
          1.96*rslt$est[[idx[k]]]$se[1,], type = "l", col = "red", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[idx[k]]]$est[3,] -
          1.96*rslt$est[[idx[k]]]$se[2,], type = "l", col = "blue", lwd = 2, lty = 2)
  lines(a.vals, rslt$est[[idx[k]]]$est[3,] + 
          1.96*rslt$est[[idx[k]]]$se[2,], type = "l", col = "blue", lwd = 2, lty = 2)
  
  if (k == 4){
    
    legend(x = 8, y = 0.02, legend=c("True ERF", "Single Imputation", "Bayesian Approach", "95% Confidence Interval"),
           col=c("darkgreen", "red", "blue", "black"),
           lty = c(1,1,1,2), lwd=2, cex=0.8)
    
  }
  
}

dev.off()

# Summary Table

tbl <- matrix(NA, nrow = length(rslt$est), ncol = 6)

for (k in 1:length(rslt$est)){
  
  bias <- round(colMeans(t(rslt$est[[k]]$bias)/rslt$est[[k]]$est[1,]), 3)
  se <- round(rowMeans(rslt$est[[k]]$se/rslt$est[[k]]$sd), 3)
  ci <- round(rowMeans(rslt$est[[k]]$cp), 3)
  
  tbl[k,] <- c(bias, sd, ci)
  
}

colnames(tbl) <- outer(names(bias), c("Bias", "SE", "CI"), FUN = "paste")[1:6]
final <- cbind(rslt$scen_idx, tbl)

write.csv(final, file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/table_2.csv")
