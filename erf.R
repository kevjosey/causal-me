# wrapper function to fit an ERF with measurement error using LOESS regression on a nonparametric models
erf <- function(a, y, x, family = gaussian(), offset = NULL, df = 4,
                a.vals = seq(min(a), max(a), length.out = 100),
                n.iter = 10000, n.adapt = 1000, thin = 10) {	
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  weights <- family$variance(family$linkinv(offset))
  n <- length(a)
  
  wrap <- bart_est(y = y, a = a, x = x, a.vals = a.vals, 
                   family = family, offset = offset,
                   n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  int.mat <- wrap$int.mat
  ybar <- family$linkinv(family$linkfun(y) - offset)
  
  # pseudo outcome
  psi <- c(ybar - muhat + mhat)
  psi[psi < 0] <- 0
  
  out <- sapply(a.vals, loess_est, psi = psi, a = a, span = 0.2, a.vals = a.vals,
                offset = offset, family = family, se.fit = TRUE, int.mat = int.mat)
  
  estimate <- out[1,]
  variance <- out[2,]
  
  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance)	
  
  return(out)
  
}

# estimate bart outcome model
bart_est <- function(a, y, x, a.vals, offset = NULL, family = gaussian(),
                     n.iter = 1000, n.adapt = 1000, thin = 10) {
  
  if (is.null(offset))
    offset <- rep(0, nrow(x))
  
  weights <- family$variance(family$linkinv(offset))
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a = a)
  colnames(xa) <- c(colnames(x), "a")
  ybar <- family$linkinv(family$linkfun(y) - offset)
  
  # for accurate simulations
  mumod <- dbarts::bart(y.train = ybar, x.train = xa, weights = weights, keeptrees = TRUE,
                        ndpost = n.iter, nskip = n.adapt, keepevery = thin, verbose = FALSE)
  muhat <- mumod$yhat.train.mean
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    # for simulations
    xa.tmp <- data.frame(x = x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa)
    return(colMeans(predict(mumod, newdata = xa.tmp, type = "ev")))
    
  })
  
  mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = weights)
  mhat <- predict(smooth.spline(a.vals, mhat.vals), x = a)$y
  
  # exposure model for integtion
  a.std <- c(c(a, a.vals) - mean(a)) / sd(a)
  dens <- density(a.std[1:n])
  phat.vals <- approx(x = dens$x, y = dens$y, xout = a.std[-(1:n)])$y / sd(a)
  
  # integration matrix
  phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)
  mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
  int.mat <- (muhat.mat - mhat.mat) * phat.mat
  
  out <- list(muhat = muhat, mhat = mhat, int.mat = int.mat)
  
  return(out)
  
}

# LOESS function
loess_est <- function(a.new, a, psi, family = gaussian(), span, offset = NULL, 
                      se.fit = FALSE, int.mat = NULL, a.vals = NULL) {
  
  if(is.null(offset))
    offset <- rep(0, times = length(a))
  
  psi <- family$linkinv(family$linkfun(psi) + offset)
  
  # construct kernel weight
  a.std <- (a - a.new)/span
  k.std <- dnorm(a.std)/span
  g.std <- cbind(1, a.std)
  
  beta <- glm(psi ~ -1 + g.std, offset = offset, weights = k.std, family = family)$coefficients
  
  # mean estimate
  mu <- family$linkinv(c(beta[1]))
  
  if (se.fit & !is.null(int.mat)) {
    
    kern.mat <- matrix(rep(dnorm((a.vals - a.new)/span)/span, n), byrow = T, nrow = n)
    g.vals <- matrix(rep(c(a.vals - a.new)/span, n), byrow = T, nrow = n)
    
    int1 <-  colMeans(k.std * t(int.mat))
    int2 <-  colMeans(a.std * k.std * t(int.mat))
    
    eta <- c(g.std %*% beta)
    w <- c(family$linkinv(eta + offset))
    z <- eta + (1/w)*(psi - family$linkinv(eta + offset))
    
    U <- solve(t(g.std) %*% diag(w * k.std) %*% g.std)
    V <- cbind(w * k.std * (z - eta) + exp(offset) * int1,
               w * a.std * k.std * (z - eta) + exp(offset) * int2)
    sig <- t(U) %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig = sig[1,1]/n))
    
  } else
    return(mu)
  
}

# k-fold cross validation to select span
cv_span <- function(a, psi, family = gaussian(), offset = NULL, 
                    folds = 5, span.seq = seq(0.05, 1, by = 0.05)) {
  
  if(is.null(offset))
    offset <- rep(0, times = length(a))
  
  n <- length(a)
  fdx <- sample(x = folds, size = min(n, 1000), replace = TRUE)
  idx <- sample(x = n, size = min(n, 1000), replace = FALSE) # for big data
  a.sub <- a[idx]
  psi.sub <- psi[idx]
  
  cv.mat <- sapply(span.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a.sub[fdx == k], loess_est, psi = psi.sub[fdx != k], a = a.sub[fdx != k], 
                      offset = offset[fdx != k], span = h, family = family, se.fit = FALSE)
      cv.vec[k] <- mean((psi.sub[fdx == k] - preds)^2, na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  span <- span.seq[which.min(cv.err)]
  
  return(span)
  
}

