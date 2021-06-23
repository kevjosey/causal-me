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
library(dbarts)
library(scales)

# Code for generating and fitting data
source("~/Github/causal-me/mclapply-hack.R")
source("~/Github/causal-me/gen-data.R")
source("~/Github/causal-me/bayes-erc.R")
source("~/Github/causal-me/bart-erc.R")
source("~/Github/causal-me/blp.R")
source("~/Github/causal-me/erc.R")

simulate <- function(scenario, n.sim, a.vals, sl.lib){
  
  # simulation arguments
  sig_gps <- scenario$sig_gps
  sig_agg <- sqrt(2)
  sig_pred <- sqrt(0.5)
  gps_scen <- scenario$gps_scen
  out_scen <- scenario$out_scen
  pred_scen <- "a"
  prob <- 0.1
  
  # gen data arguments
  n <- scenario$n # c(500, 800)
  mult <- scenario$mult # c(100, 200)
  
  # gibbs sampler stuff
  thin <- 10
  n.iter <- 2000
  n.adapt <- 500
  h.a <- 0.5
  h.gamma <- 0.03
  scale <- 1e6
  shape <- rate <- 1e-3
  deg.num <- 2
  
  # dr arguments
  span <- ifelse(n == 800, 0.125, 0.25)
  family <- poisson()
  
  # initialize output
  est <- array(NA, dim = c(n.sim, 6, length(a.vals)))
  se <- array(NA, dim = c(n.sim, 5, length(a.vals)))
  
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
    z_hat <- aggregate(s_hat, by = list(s.id), mean)[,2]
    
    # real
    obs_hat <- try(erc(y = y, a = a, x = x, offset = offset, family = family,
                       a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    
    # naive
    naive_hat <- try(erc(y = y, a = z_hat, x = x, offset = offset, family = family,
                     a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num))
    
    # blp
    blp_hat <- try(erc(y = y, a = a_hat, x = x, offset = offset, family = family,
                     a.vals = a.vals, sl.lib = sl.lib, span = span, deg.num = deg.num), silent = TRUE)
    
    # Bayesian Approach
    bayes_hat <- try(bayes_erc(s = s, star = s_tilde, y = y, offset = offset,
                            s.id = s.id, id = id, w = w, x = x, family = family,
                            a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                            n.iter = n.iter, n.adapt = n.adapt, thin = thin, 
                            h.a = h.a, h.gamma = h.gamma, deg.num = deg.num), silent = TRUE)
    
    # BART Approach
    bart_hat <- try(bart_erc(s = s, star = s_tilde, y = y, offset = offset,
                             s.id = s.id, id = id, w = w, x = x, family = family, 
                             a.vals = a.vals, span = span, scale = scale, shape = shape, rate = rate,
                             h.a = h.a, n.iter = n.iter, n.adapt = n.adapt, thin = thin), silent = TRUE)
    
    # estimates
    est <- rbind(predict_example(a = a.vals, x = x, out_scen = out_scen),
                 if (!inherits(obs_hat, "try-error")) {obs_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(naive_hat, "try-error")) {naive_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(blp_hat, "try-error")) {blp_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bayes_hat, "try-error")) {bayes_hat$estimate} else {rep(NA, length(a.vals))},
                 if (!inherits(bart_hat, "try-error")) {bart_hat$estimate} else {rep(NA, length(a.vals))})
    
    #standard error
    se <- rbind(if (!inherits(obs_hat, "try-error")) {sqrt(obs_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(naive_hat, "try-error")) {sqrt(bayes_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(blp_hat, "try-error")) {sqrt(blp_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bayes_hat, "try-error")) {sqrt(bayes_hat$variance)} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {sqrt(bart_hat$variance)} else {rep(NA, length(a.vals))})
    
    # coverage probability
    cp <- rbind(if (!inherits(obs_hat, "try-error")) {as.numeric((est[2,] - 1.96*se[1,]) < est[1,] & (est[2,] + 1.96*se[1,]) > est[1,])} else {rep(NA, length(a.vals))},
                if (!inherits(naive_hat, "try-error")) {as.numeric((est[3,] - 1.96*se[2,]) < est[1,] & (est[3,] + 1.96*se[2,]) > est[1,])} else {rep(NA, length(a.vals))},
                if (!inherits(blp_hat, "try-error")) {as.numeric((est[4,] - 1.96*se[3,]) < est[1,] & (est[4,] + 1.96*se[3,]) > est[1,])} else {rep(NA, length(a.vals))},
                if (!inherits(bayes_hat, "try-error")) {as.numeric((est[5,] - 1.96*se[4,]) < est[1,] & (est[5,] + 1.96*se[4,]) > est[1,])} else {rep(NA, length(a.vals))},
                if (!inherits(bart_hat, "try-error")) {as.numeric((est[6,] - 1.96*se[5,]) < est[1,] & (est[6,] + 1.96*se[5,]) > est[1,])} else {rep(NA, length(a.vals))})
    
    return(list(est = est, se = se, cp = cp))
     
  }, mc.cores = 25)
  
  est <- abind(lapply(out, function(lst, ...) lst$est), along = 3)
  se <- abind(lapply(out, function(lst, ...) lst$se), along = 3)
  cp <- abind(lapply(out, function(lst, ...) lst$cp), along = 3)
  
  out_est <- t(apply(est, 1, rowMeans, na.rm = T))
  colnames(out_est) <- a.vals
  rownames(out_est) <- c("ERF","DR","Naive","BLP","Bayes","BART")
  
  compare <- matrix(est[1,,], nrow = length(a.vals), ncol = n.sim)
  out_bias <- t(apply(est[2:6,,], 1, function(x) rowMeans(abs(x - compare), na.rm = T)))
  colnames(out_bias) <- a.vals
  rownames(out_bias) <- c("DR","Naive","BLP","Bayes","BART")
  
  out_sd <- t(apply(est[2:6,,], 1, function(x) apply(x, 1, sd, na.rm = T)))
  colnames(out_sd) <- a.vals
  rownames(out_sd) <- c("DR","Naive","BLP","Bayes","BART")
  
  out_se <- t(apply(se, 1, rowMeans, na.rm = T))
  colnames(out_se) <- a.vals
  rownames(out_se) <- c("DR","Naive","BLP","Bayes","BART")
  
  out_cp <- t(apply(cp, 1, rowMeans, na.rm = T))
  colnames(out_cp) <- a.vals
  rownames(out_se) <- c("DR","Naive","BLP","Bayes","BART")
  
  return(list(est = out_est, bias = out_bias, sd = out_sd, se = out_se, cp = out_cp))
  
}

# for replication
set.seed(42)

# simulation scenarios
a.vals <- seq(6, 10, by = 0.04)
sl.lib <- c("SL.mean","SL.glm")
n.sim <- 1000

n <- c(400, 800)
mult <- c(5, 10)
gps_scen <- c("a", "b")
out_scen <- c("a", "b")
sig_gps <- c(1, sqrt(4))

scen_mat <- expand.grid(n = n, mult = mult, gps_scen = gps_scen, out_scen = out_scen, sig_gps = sig_gps)
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
pdf(file = filename, width = 9, height = 9)
par(mfrow = c(2,2))

for (k in 1:4){
  
  plot(a.vals, rslt$est[[idx[k]]]$est[1,], type = "l", col = hue_pal()(6)[1], lwd = 2,
       xlab = "Exposure", ylab = "Rate of Event", main = plotnames[k],
       ylim = c(0,0.12))
  grid(lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[2,], type = "l", col = hue_pal()(6)[2], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[3,], type = "l", col = hue_pal()(6)[3], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[4,], type = "l", col = hue_pal()(6)[4], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[5,], type = "l", col = hue_pal()(6)[5], lwd = 2, lty = 1)
  lines(a.vals, rslt$est[[idx[k]]]$est[5,], type = "l", col = hue_pal()(6)[6], lwd = 2, lty = 1)
  
  if (k == 4){
    
    legend(6, 0.15, legend=c("True ERF", "Observed", "Naive", "BLP", "BART"),
           col=hue_pal()(5),
           lty = c(1,1,1,1,1), lwd=2, cex=0.8)
    
    
  }
  
}

dev.off()

# Summary Table

tbl <- matrix(NA, nrow = length(rslt$est), ncol = 6)

for (k in 1:length(rslt$est)){
  
  bias <- round(colMeans(t(rslt$est[[k]]$bias)/rslt$est[[k]]$est[1,]), 3)
  se <- round(rowMeans(rslt$est[[k]]$se/rslt$est[[k]]$sd), 3)
  ci <- round(rowMeans(rslt$est[[k]]$cp), 3)
  
  tbl[k,] <- c(bias, se, ci)
  
}

colnames(tbl) <- outer(names(bias), c("Bias", "SE", "CI"), FUN = "paste")[1:6]
final <- cbind(rslt$scen_idx, tbl)

write.csv(final, file = "~/Dropbox (Personal)/Projects/ERC-EPE/Output/table_2.csv")
