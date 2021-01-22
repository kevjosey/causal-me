library(nnet)
library(MASS)
library(mvtnorm)
library(npcausal)

source("D:/Github/causal-me/rc-sim/fun.R")
source("D:/Github/causal-me/mclapply-hack.R")

scenarios <- 1
n.sim <- 100

sim.fun<-function(scenarios = 1, n.sim = 100, ...){
  
  #gamma1<-NULL
  ########default
  if (scenarios==1){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1"}
  ########correlation low
  if (scenarios==2){text<-"phi=0.2,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1"}
  ######## poor model fit
  if (scenarios==3){text<-"phi=0.8,sd_rc=10,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=1"}
  ######## outcome effect small
  if (scenarios==4){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=0.5,beta_par=1,correct_RC=1"}
  ####### confounder on outcome large
  if (scenarios==5){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=5,correct_RC=1"}
  ####### confounder on exposure large
  if (scenarios==6){text<-"phi=0.8,sd_rc=1,tau_par=1.6,beta_1=1,beta_par=1,correct_RC=1"}
  ####### RC model misspecified
  if (scenarios==7){text<-"phi=0.8,sd_rc=1,tau_par=0.8,beta_1=1,beta_par=1,correct_RC=0"}
  
  est <- array(NA, dim = c(n.sim, 4, length(a.vals)))
  se <- array(NA, dim = c(n.sim, 3, length(a.vals)))
  degree <- 3
  a.rng <- c(-20, 20)
  bw.seq <- seq(0.2, 5, by = 0.1)
  n.boot <- 20
  
  for (i in 1:n.sim){
    
    simulated.data<-eval(parse(text=paste0("data.generate(sample_size=2000,seed = i,",text,")")))
    validation<-simulated.data[sample.int(nrow(simulated.data),500),]
    
    y <- simulated.data$Y
    a <- simulated.data$treat
    x <- simulated.data[,c("cf1", "cf2", "cf3", "cf4", "cf5", "cf6")]
    
    rc.naive <- lm(treat ~ treat.w, validation)
    z.naive <- predict(rc.naive, newdata = simulated.data)
    
    rc.model <- lm(treat~treat.w+cova1+cova2+cova3, validation)
    z.rc <- predict(rc.model, newdata = simulated.data)
    sigma <- sigma(rc.model)
    
    naive <- ctseff(a = z.naive, y = y, x = x, a.rng = a.rng, n.pts = n.pts, bw.seq = bw.seq, sl.lib = sl.lib)
    
    naive_est <- naive$res$est
    naive_var <- naive$res$se^2
    
    rc <- ctseff(a = z.rc, y = y, x = x, a.rng = a.rng, n.pts = n.pts, bw.seq = bw.seq, sl.lib = sl.lib)
    
    rc_est <- rc$res$est
    rc_var <- rc$res$se^2
    
    simex <- simex_cts(z = z.rc, y = y, x = x, a.rng = a.rng, bw.seq = bw.seq,
                       n.pts = n.pts, sigma = sigma, n.boot = n.boot, degree = degree)
    
    simex_est <- simex$estimate
    simex_var <- simex$variance
    
    n <- dim(x)[1]
    a.min <- a.rng[1]
    a.max <- a.rng[2]
    a.vals <- seq(a.min, a.max, length.out = n.pts)
    
    beta<-c(1,3,2,1,4,2,1,-1,1,-1)
    if (scenarios == 4)
      beta<-c(0.5,3,2,1,4,2,1,-1,1,-1)
    if (scenarios == 5)
      beta<-c(1,c(3,2,1,4,2,1,-1,1,-1)*5)
    #produce outcome Y
    mu <- sapply(a.vals, function(z, ...){
      mean(as.matrix(cbind(z, x, z*x[,c(2,4,6)]))%*%beta)
    })
    
    # combine output
    est[i,1,] <- mu
    est[i,2:4,] <- rbind(naive_est, rc_est, simex_est)
    se[i,1:3,] <- sqrt(rbind(naive_var, rc_var, simex_var))

  }

  return(list(est = est, se))
  
}

out <- mclapply.hack(1:7,sim.fun, n.sim = n.sim, mc.cores=7)
