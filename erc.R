
# wrapper function to fit a hierarchical, doubly-robust ERC using LOESS regression on a nonparametric model
erc <- function(a, y, x, family = gaussian(), offset = rep(0, length(a)),
                a.vals = seq(min(a), max(a), length.out = 100), deg.num = 2,
                span = NULL, span.seq = seq(0.05, 1, by = 0.05), k = 5,
                sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.earth")){	
  
  n <- length(a)
  
  wrap <- np_est(y = y, a = a, x = x, offset = offset, deg.num = deg.num,
                 a.vals = a.vals, family = family, sl.lib = sl.lib)
  
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  pihat <- wrap$pihat
  phat <- wrap$phat
  int <- wrap$int
  y_ <- family$linkinv(family$linkfun(y) - offset)
  psi <- c((y_ - muhat) + mhat)
  
  if(is.null(span)) {
    
    idx <- sample(x = n, size = min(n, 1000), replace = FALSE)
    
    a.sub <- a[idx]
    psi.sub <- psi[idx]
    int.sub <- int[idx]
    
    folds <- sample(x = k, size = min(n, 1000), replace = TRUE)
    
    cv.mat <- sapply(span.seq, function(h, ...) {
      
      cv.vec <- rep(NA, k)
      
      for(j in 1:k) {
        
        preds <- sapply(a.sub[folds == j], dr_est, psi = psi.sub[folds != j], a = a.sub[folds != j], 
                        int = int.sub[folds != j], span = h, family = family, se.fit = FALSE)
        cv.vec[j] <- mean((psi.sub[folds == j] - preds)^2, na.rm = TRUE)
        
      }
      
      return(cv.vec)
      
    })
    
    cv.err <- colMeans(cv.mat)
    span <- span.seq[which.min(cv.err)]
    
  }
  
  dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = int, 
                   span = span, family = gaussian(), se.fit = TRUE)
  
  estimate <- dr_out[1,]
  variance <- dr_out[2,]
  
  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance, span = span)	
  
  return(out)
  
}

# LOESS function
dr_est <- function(newa, a, psi, int, span, family = gaussian(), se.fit = FALSE) {
  
  a.std <- a - newa
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  a.std <- a.std[idx]
  psi <- psi[idx]
  max.a.std <- max(abs(a.std))
  k.std <- c((1 - abs(a.std/max.a.std)^3)^3)
  gh <- cbind(1, a.std)
  b <- optim(par = c(0,0), fn = opt_fun, k.std = k.std, psi = psi, gh = gh, family = family)
  mu <- family$linkinv(c(b$par[1]))
  
  if (se.fit){
    
    int <- int[idx]
    gh.inv <- solve(t(gh) %*% diag(k.std) %*% gh)
    v.inf <- (psi + int - family$linkinv(c(gh%*%b$par)))^2
    sig <- gh.inv %*% t(gh) %*% diag(k.std) %*% diag(v.inf) %*% diag(k.std) %*% gh %*% gh.inv
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)
  
}

np_est <- function(a, y, x, a.vals = a.vals, family = gaussian(), offset = rep(0, length(a)), deg.num = 2,
                   sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.earth")) {
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a = a)
  colnames(xa) <- c(colnames(x), "a")
  y_ <- family$linkinv(family$linkfun(y) - offset)
  
  # estimate nuisance GPS functions via super learner
  pimod <- SuperLearner(Y = a, X = x, family = gaussian(), SL.library = sl.lib)
  pimod.vals <- c(pimod$SL.predict)
  pi2mod <- SuperLearner(Y = (a - pimod.vals)^2, X = x, family = gaussian(), SL.library = sl.lib)
  pi2mod.vals <- c(pi2mod$SL.predict)
  pi2mod.vals[pi2mod.vals <= 0] <- .Machine$double.eps
  
  # exposure models
  pihat <- dnorm(a, pimod.vals, sqrt(pi2mod.vals))
  phat.vals <- sapply(a.vals, function(a.tmp, ...) 
    mean(dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals))))
  phat <- predict(smooth.spline(a.vals, phat.vals), x = a)$y
  phat[which(phat < 0)] <- .Machine$double.eps
  phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)
  
  # for accurate simulations
  mumod <- glm(y ~ . - a + poly(a, deg.num) + a:x1, data = xa, family = family, offset = offset)
  muhat <- predict(mumod, newdata = xa, type = "response")
  
  # estimate nuisance outcome model with SuperLearner
  # mumod <- SuperLearner(Y = y_, X = xa, family = gaussian(), SL.library = sl.lib)
  # muhat <- c(mumod$SL.predict)
  
  # predict marginal outcomes given a.vals (or a.agg)
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    # for simulations
    xa.tmp <- data.frame(x = x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa)
    return(predict(mumod, newdata = xa.tmp, type = "response"))
    
    # general approach corresponding to SuperLearner model
    # xa.tmp <- data.frame(x = x, a = a.tmp)
    # colnames(xa.tmp) <- colnames(xa)
    # return(predict(mumod, newdata = xa.tmp)$pred)
    
  })
  
  # aggregate muhat.vals and integrate for influence function
  mhat.vals <- colMeans(muhat.mat)
  mhat <- predict(smooth.spline(a.vals, mhat.vals), x = a)$y
  mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
  
  # integrate
  intfn <- (muhat.mat - mhat.mat) * phat.mat
  int <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
                 (intfn[,-1] + intfn[,-length(a.vals)]) / 2, 1, sum)
  
  out <- list(muhat = muhat, mhat = mhat, pihat = pihat, phat = phat, int = int, 
              mhat.vals = mhat.vals, phat.vals = phat.vals)
  
  return(out)
  
}

# optimization used in dr_est
opt_fun <- function(par, k.std, psi, gh, family) {
  
  sum(k.std*(psi - family$linkinv(c(gh %*% par)))^2)
  
}

