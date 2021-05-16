ipw <- function(a, x, beta, sigma2, a.vals) {
  
  n <- length(a)
  x.new <- rbind(x, x[rep(1:n, length(a.vals)), ])
  a.new <- c(a, rep(a.vals, each = n))
  pimod.vals <- c(x.new %*% beta)
  pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2))
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
  phat[which(phat < 0)] <- 1e-6
  out <- phat/pihat
  return(out)
  
}

marg <- function(a, x, nsa, gamma, a.vals) {
  
  xa <- cbind(x, predict(nsa, a))
  
  xa.new.list <- lapply(a.vals, function(a.tmp, ...) {
    
    cbind(x, matrix(rep(c(predict(nsa, a.tmp)), n), byrow = TRUE, nrow = n))
    
  })
  
  xa.new <- rbind(xa, do.call(rbind, xa.new.list))
  colnames(xa.new) <- colnames(xa)
  muhat.vals <- family$linkinv(c(xa.new %*% gamma))
  muhat <- muhat.vals[1:n]
  muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
  return(list(muhat = muhat, mhat = mhat))
  
} 

split.along.dim <- function(a, n) {
  
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
  
}