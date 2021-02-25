#function for generating simulated data, including both "true" and "error-prone" exposurs
data.generate<-function(sample_size=2000,seed=300,sd=20,phi=0.8,sd_rc=1,tau_par=0.6,beta_1=1,beta_par=1,correct_RC=1){
  
  options(digits=4) # only print 4 sig digits
  set.seed(seed)
  size<-sample_size
  #pre-treatment variables (confounders)
  cf1<-mvrnorm(n=size,mu=c(0,0,0),Sigma=matrix(c(2,1,-1,1,1,-0.5,-1,-0.5,1),ncol=3))
  cf2<-sample(c((-2):2),size,replace = T)
  cf3<-runif(size,min=-3,max=3)
  cf4<-rchisq(size,1)
  cf<-cbind(cf0=1,cf1,cf2,cf3,cf4)
  #coefficients of linear regression model
  tau<-(tau_par*c(1,1,2,1.5,3,2,3))
  
  # covariates in the RC model
  cova1<-cf1[,1]
  cova2<-rnorm(size,0,sd=4)
  cova3<-runif(size,min=-5,max=5)
  cova<-cbind(cova1,cova2,cova3)
  zeta<-c(2,1,3,phi) 
  
  treat.w<-as.numeric()
  treat<-as.numeric()
  for (i in 1:size){
    treat[i]<-sum(cf[i,]*tau)+rnorm(1,mean=0,sd=sd)
    if (correct_RC==1){
      treat.w[i]<-sum(c(cova[i,],treat[i])*zeta) + rnorm(1,mean=0,sd=sd_rc) 
    } else if (correct_RC==0){
      treat.w[i]<-sum(c(cova[i,],treat.w[i])*zeta)-
        5*log(abs(cova2[i]))-3*cova1[i]^2-20*sin(cova3[i])+
        5*exp(rnorm(1,0,1))+rnorm(1,mean=0,sd=sd_rc) 
    }
  }
  
  #coefficients of linear model
  beta<-c(beta_1,c(3,2,1,4,2,1,-1,1,-1)*beta_par)
  #produce outcome Y
  Y<-as.numeric()
  for (i in 1:size){
    Y[i]<-sum(c(treat[i],cf[i,-1], treat[i]*cf[i,c(2,4,6)])*beta) + rnorm(1, mean = 0, sd = 1)
  }
  
  #treat<-treat
  #treat.w<-treat.w+100
  simulated.data<-data.frame(cbind(Y,treat,cf[,-1],cova,treat.w))
  colnames(simulated.data)[3:8]<-c("cf1","cf2","cf3","cf4","cf5","cf6")
  return(simulated.data)
  
}

# needed to add an na.rm = TRUE
ctseff <- function (y, a, x, bw.seq, n.pts = 100, a.rng = c(min(a), max(a)), 
                    sl.lib = c("SL.earth", "SL.glm", "SL.glm.interaction")) {
  
  require("SuperLearner")
  require(KernSmooth)
  
  kern <- function(t) {
    dnorm(t)
  }
  
  n <- dim(x)[1]
  a.min <- a.rng[1]
  a.max <- a.rng[2]
  a.vals <- seq(a.min, a.max, length.out = n.pts)
  xa.new <- rbind(cbind(x, a),
                  cbind(x[rep(1:n, length(a.vals)),],
                        a = rep(a.vals, rep(n, length(a.vals)))))
  x.new <- xa.new[, -dim(xa.new)[2]]
  x <- data.frame(x)
  x.new <- data.frame(x.new)
  colnames(x) <- colnames(x.new)
  xa.new <- data.frame(xa.new)
  
  pimod <- SuperLearner(Y = a, X = data.frame(x), SL.library = sl.lib, newX = x.new)
  pimod.vals <- pimod$SL.predict
  pi2mod <- SuperLearner(Y = (a - pimod.vals[1:n])^2, X = x,  SL.library = sl.lib, newX = x.new)
  pi2mod.vals <- pi2mod$SL.predict
  
  mumod <- SuperLearner(Y = y, X = cbind(x, a), SL.library = sl.lib, newX = xa.new)
  muhat.vals <- mumod$SL.predict
  
  a.std <- (xa.new$a - pimod.vals)/sqrt(pi2mod.vals)
  pihat.vals <- approx(density(a.std)$x, density(a.std[1:n])$y, xout = a.std)$y/sqrt(pi2mod.vals)
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  varpihat <- predict(smooth.spline(a.vals, apply(pihat.mat, 2, mean)), x = a)$y
  varpihat.mat <- matrix(rep(apply(pihat.mat, 2, mean), n), byrow = T, nrow = n)
  
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  mhat <- predict(smooth.spline(a.vals, apply(muhat.mat, 2, mean)), x = a)$y
  mhat.mat <- matrix(rep(apply(muhat.mat, 2, mean), n), byrow = T, nrow = n)
  
  pseudo.out <- (y - muhat)/(pihat/varpihat) + mhat
  
  w.fn <- function(bw) {
    w.avals <- NULL
    for (a.val in a.vals) {
      a.std <- (a - a.val)/bw
      kern.std <- kern(a.std)/bw
      w.avals <- c(w.avals, mean(a.std^2 * kern.std) * (kern(0)/bw)/
                     (mean(kern.std) * mean(a.std^2 * kern.std) - mean(a.std * kern.std)^2))
    }
    return(w.avals/n)
  }
  
  hatvals <- function(bw) {
    approx(a.vals, w.fn(bw), xout = a)$y
  }
  
  cts.eff.fn <- function(out, bw) {
    approx(locpoly(a, out, bandwidth = bw), xout = a)$y
  }
  
  risk.fn <- function(h) {
    hats <- hatvals(h)
    mean(((pseudo.out - cts.eff.fn(pseudo.out, bw = h))/(1 - hats))^2, na.rm = TRUE)
  }
  
  risk.est <- sapply(bw.seq, risk.fn)
  h.opt <- bw.seq[which.min(round(risk.est, 3))]
  bw.risk <- data.frame(bw = bw.seq, risk = risk.est)
  est <- approx(locpoly(a, pseudo.out, bandwidth = h.opt), xout = a.vals)$y
  se <- NULL
  
  for (a.val in a.vals) {
    a.std <- (a - a.val)/h.opt
    kern.std <- kern(a.std)/h.opt
    beta <- coef(lm(pseudo.out ~ a.std, weights = kern.std))
    Dh <- matrix(c(mean(kern.std), mean(kern.std * a.std), 
                   mean(kern.std * a.std), mean(kern.std * a.std^2)), nrow = 2)
    kern.mat <- matrix(rep(kern((a.vals - a.val)/h.opt)/h.opt, n), byrow = T, nrow = n)
    g2 <- matrix(rep((a.vals - a.val)/h.opt, n), byrow = T, 
                 nrow = n)
    intfn1.mat <- kern.mat * (muhat.mat - mhat.mat) * varpihat.mat
    intfn2.mat <- g2 * kern.mat * (muhat.mat - mhat.mat) * 
      varpihat.mat
    int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n),
                         byrow = T, nrow = n) * 
                    (intfn1.mat[, -1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
    int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n),
                         byrow = T, nrow = n) * 
                    (intfn2.mat[, -1] + intfn2.mat[, -length(a.vals)])/2, 1, sum)
    sigma <- cov(t(solve(Dh) %*% rbind(kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int1,
                                       a.std * kern.std * (pseudo.out - beta[1] - beta[2] * a.std) + int2)))
    se <- c(se, sqrt(sigma[1, 1]))
  }
  
  ci.ll <- est - 1.96 * se/sqrt(n)
  ci.ul <- est + 1.96 * se/sqrt(n)
  res <- data.frame(a.vals, est, se, ci.ll, ci.ul)
  return(invisible(list(res = res, bw.risk = bw.risk)))
  
}

simex_cts <- function(z, y, x, sigma, n.boot = 100, degree = 2, lambda = seq(0.1, 2.1, by = 0.25),
                      a.rng = seq(min(a), max(a)), bw.seq = seq(.2, 2, length.out = 100), n.pts = 100,
                      sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction")) {
  
  eps <- replicate(n.boot, rnorm(length(z), 0, sigma))
  
  l.vals <- mclapply.hack(lambda, function(lam, ...){
    
    z.mat <- z + sqrt(lam)*eps
    
    vals <- apply(z.mat, 2, ctseff, y = y, x = x, a.rng = a.rng, 
                  n.pts = n.pts, bw.seq = bw.seq, sl.lib = sl.lib)
    
    mu.vals <- matrix(unlist(lapply(vals, function(x) x$res$est)), ncol = length(vals))
    sig.vals <- matrix(unlist(lapply(vals, function(x) x$res$se^2)), ncol = length(vals))
    
    m.vals <- matrix(rep(rowMeans(mu.vals), n.boot), nrow = length(a.vals), ncol = n.boot)
    s.hat <- rowSums((mu.vals - m.vals)^2)/(n.boot - 1)
    tau.hat <- rowMeans(sig.vals)
    
    return(list(estimate = rowMeans(mu.vals), variance = tau.hat - s.hat))
    
  })
  
  if (any(lambda == 0)){
    
    Psi <- t(matrix(unlist(lapply(l.vals, function(x) x$estimate)), ncol = length(l.vals)))[,-which(lambda == 0)]
    Phi <- t(matrix(unlist(lapply(l.vals, function(x) x$variance)), ncol = length(l.vals)))[,-which(lambda == 0)]
    
  } else {
    
    Psi <- t(matrix(unlist(lapply(l.vals, function(x) x$estimate)), ncol = length(l.vals)))
    Phi <- t(matrix(unlist(lapply(l.vals, function(x) x$variance)), ncol = length(l.vals)))
    
  }
  
  L <- cbind(1, poly(lambda, degree = degree, raw = TRUE))
  chi <- c(1, poly(-1, degree = degree, raw = TRUE))
  estimate <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Psi)
  variance <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Phi)
  
  return(list(estimate = estimate, variance = variance, Psi = Psi, Phi = Phi, lambda = lambda))
  
}
