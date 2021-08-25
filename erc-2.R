
# wrapper function to fit a hierarchical, doubly-robust ERC using LOESS regression on a nonparametric model
erc <- function(a, y, x, family = gaussian(), offset = NULL, weights = NULL,
                a.vals = seq(min(a), max(a), length.out = 100), deg.num = 2,
                span = 0.5, span.seq = seq(0.05, 1, by = 0.05), k = 5){	
  
  
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  n <- length(a)
  
  wrap <- glm_est(y = y, a = a, x = x, offset = offset, deg.num = deg.num,
                 a.vals = a.vals, family = family)
  
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  int.mat <- wrap$int.mat
  
  y_ <- family$linkinv(family$linkfun(y) - offset)
  psi <- c((y_ - muhat) + mhat)
  
  # if(is.null(span)) {
  #   
  #   idx <- sample(x = n, size = min(n, 1000), replace = FALSE)
  #   
  #   a.sub <- a[idx]
  #   psi.sub <- psi[idx]
  #   int.sub <- int[idx]
  #   
  #   folds <- sample(x = k, size = min(n, 1000), replace = TRUE)
  #   
  #   cv.mat <- sapply(span.seq, function(h, ...) {
  #     
  #     cv.vec <- rep(NA, k)
  #     
  #     for(j in 1:k) {
  #       
  #       preds <- sapply(j, a.sub, dr_est, psi = psi.sub[folds != j], a = a.sub[folds != j], 
  #                       int = int.sub[folds != j], span = h, family = family, se.fit = FALSE)
  #       cv.vec[j] <- mean((psi.sub[folds == j] - preds)^2, na.rm = TRUE)
  #       
  #     }
  #     
  #     return(cv.vec)
  #     
  #   })
  #   
  #   cv.err <- colMeans(cv.mat)
  #   span <- span.seq[which.min(cv.err)]
  #   
  # }
  
  dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, family = gaussian(), 
                   span = span, int.mat = int.mat, se.fit = TRUE)
  
  estimate <- dr_out[1,]
  variance <- dr_out[2,]
  
  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance, span = span)	
  
  return(out)
  
}

# LOESS function
dr_est <- function(newa, a, psi, span, family = gaussian(), se.fit = FALSE, int.mat = NULL) {
  
  a.std <- a - newa
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  a.std <- a.std[idx]
  psi <- psi[idx]
  max.a.std <- max(abs(a.std))
  k.std <- c((1 - abs(a.std/max.a.std)^3)^3)
  gh <- cbind(1, a.std)
  bh <- optim(par = c(0,0), fn = opt_fun, k.std = k.std, psi = psi, gh = gh, family = family)
  mu <- family$linkinv(c(bh$par[1]))
  
  if (se.fit & !is.null(int.mat)){
    
    Dh <- matrix(c(mean(k.std), mean(k.std * a.std),
                   mean(k.std * a.std), mean(k.std * a.std^2)), nrow = 2)
    
    intfn1.mat <- k.std*int.mat[idx,]
    intfn2.mat <- gh[,2]*k.std*int.mat[idx,]
    int1 <- rowMeans(intfn1.mat)
    int2 <- rowMeans(intfn2.mat)
    
    sigma <- cov(t(solve(Dh) %*% rbind(k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int1,
                                       a.std * k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int2)))
    
    return(c(mu = mu, sig = sigma[1,1]))
    
  } else
    return(mu)
  
}

glm_est <- function(a, y, x, a.vals, family = gaussian(), offset = rep(0, length(a)), deg.num = 2) {
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a = a)
  colnames(xa) <- c(colnames(x), "a")
  y_ <- family$linkinv(family$linkfun(y) - offset)
  
  # for accurate simulations
  mumod <- glm(y ~ . - a + poly(a, deg.num) + a:x1, data = xa, family = family, offset = offset)
  muhat <- predict(mumod, newdata = xa, type = "response")
  
  # muhat.mat <- sapply(a.vals, function(a.tmp, ...){
  #   xa.tmp <- data.frame(x = x, a = a.tmp)
  #   colnames(xa.tmp) <- colnames(xa)
  #   return(predict(mumod, newdata = xa.tmp, type = "response"))
  # })
  # 
  # mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
  # int.mat <- t(apply(muhat.mat, 1, function(val,...) 
  #   predict(smooth.spline(a.vals, val), x = a)$y)) - 
  #   matrix(rep(mhat, n), byrow = T, nrow = n)
  
  muhat.mat <- sapply(a, function(a.tmp, ...) {
    
    # for simulations
    xa.tmp <- data.frame(x = x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa)
    return(predict(mumod, newdata = xa.tmp, type = "response"))
    
  })
  
  mhat <- colMeans(muhat.mat)
  int.mat <- muhat.mat - matrix(rep(mhat, n), byrow = T, nrow = n)
  
  out <- list(muhat = muhat, mhat = mhat, int.mat = int.mat)
  
  return(out)
  
}

# optimization used in dr_est
opt_fun <- function(par, k.std, psi, gh, family) {
  
  sum(k.std*(psi - family$linkinv(c(gh %*% par)))^2)
  
}
