
# wrapper function to fit a hierarchical, doubly-robust ERC using LOESS regression on a nonparametric model
erc <- function(a, y, x, wts,
                a.vals = seq(min(a), max(a), length.out = 20), 
                span = NULL, span.seq = seq(0.15, 1, by = 0.05), k = 10) {	
  
  m <- length(a)
  
  wrap <- np_est(y = y, a = a, x = x, wts = wts)
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  pihat <- wrap$pihat
  phat <- wrap$phat
  int <- wrap$int
  psi <- (y - wts*muhat)/(pihat/phat) + mhat
  
  if(is.null(span)) {
    
    folds <- sample(x = k, size = m, replace = TRUE)
    
    cv.mat <- sapply(span.seq, function(h, ...) {
      
      cv.vec <- rep(NA, k)
      
      for(j in 1:k) {
        
        preds <- sapply(a[folds == j], dr_est, psi = psi[folds != j], a = a[folds != j], 
                        int = int[folds != j],  span = h, se.fit = FALSE)
        cv.vec[j] <- mean((psi[folds == j] - preds)^2, na.rm = TRUE)
        # some predictions result in `NA` because of the `x` ranges in each fold
        
      }
      
      return(cv.vec)
      
    })
    
    cv.err <- colMeans(cv.mat)
    span <- span.seq[which.min(cv.err)]
    
  }
  
  dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = int, span = span, se.fit = TRUE)
  
  estimate <- dr_out[1,]
  variance <- dr_out[2,]
  
  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance)	
  
  return(out)
  
}

# LOESS function
dr_est <- function(newa, a, psi, int, span, se.fit = FALSE) {
  
  a.std <- a - newa
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  knn <- rep(0, length(a))
  knn[idx] <- 1
  max.a.std <- max(abs(a.std*knn))
  k.std <- ((70/81)*(1 - abs(a.std/max.a.std)^3)^3)*knn
  
  # a.std <- (a - newa) / h
  # k.std <- dnorm(a.std) / h
  
  gh <- cbind(1, a.std)
  gh.inv <- solve(t(gh) %*% diag(k.std) %*% gh)
  b <- optim(par = c(0,0), fn = opt_fun, k.std = k.std, psi = psi, gh = gh)
  mu <- plogis(b$par[1])
  
  if (se.fit){
    
    v.inf <- (psi + int - plogis(c(gh%*%b$par)))^2
    sig <- gh.inv %*% t(gh) %*% diag(k.std) %*% diag(v.inf) %*% diag(k.std) %*% gh %*% gh.inv
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)
  
}

# Nonparametric estimation
np_est <- function(a, y, x, wts, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.earth")) {
  
  # set up evaluation points & matrices for predictions
  m <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(y = y, a = a, x = x)
  xa.new <- data.frame(a = a, x = x, wts = 1)
  colnames(xa.new) <- c(colnames(xa)[-1], "wts")
  
  # estimate nuisance outcome model with SuperLearner
  # mumod <- SuperLearner(Y = y, X = x, family = poisson(), SL.library = sl.lib)
  mumod <- glm(y ~ .^2, offset = log(c(wts)), family = poisson(log), data = xa)
  muhat <- predict(mumod, type = "response", newdata = xa.new)
  
  # estimate nuisance GPS functions via super learner
  pimod <- SuperLearner(Y = a, X = x, family = gaussian(), SL.library = sl.lib)
  pimod.vals <- c(pimod$SL.predict)
  pi2mod <- SuperLearner(Y = (a - pimod.vals)^2, X = x, 
                         family = gaussian(), SL.library = "SL.mean")
  pi2mod.vals <- c(pi2mod$SL.predict)
  
  # exposure models
  pihat <- dnorm(a, pimod.vals, sqrt(pi2mod.vals))
  phat <- sapply(a, function(a.tmp, ...) weighted.mean(dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals))), w = wts)
  
  # predict marginal outcomes given a.vals (or a.agg)
  muhat.mat <- sapply(a, function(a.tmp, ...) {
    
    xa.tmp <- data.frame(a = a.tmp, x = x, wts = 1)
    colnames(xa.tmp) <- c(colnames(xa)[-1], "wts") 
    return(c(predict(mumod, newdata = xa.tmp, type = "response")))
    
  })
  
  # aggregate muhat.vals and integrate for influence
  mhat <- colMeans(muhat.mat)
  mhat.mat <- matrix(rep(mhat, m), byrow = TRUE, nrow = m)
  int <- rowMeans(muhat.mat - mhat.mat)
  
  out <- list(muhat = muhat, mhat = mhat, pihat = pihat, phat = phat,int = int)
  
  return(out)
  
}

# optimization used in dr_est
opt_fun <- function(par, k.std, psi, gh) {
  
  sum(k.std*(psi - exp(c(gh %*% par)))^2)
  
}