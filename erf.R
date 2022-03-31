# wrapper function to fit an ERF with measurement error using LOESS regression on a nonparametric models
erf <- function(a, y, x, family = gaussian(), offset = NULL, weights = NULL,
                a.vals = seq(min(a), max(a), length.out = 100),
                n.iter = 10000, n.adapt = 1000, thin = 10, span = NULL, 
                span.seq = seq(0.05, 1, by = 0.05), k = 5) {	
  
  
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  n <- length(a)
  
  wrap <- bart_est(y = y, a = a, x = x, a.vals = a.vals,
                   family = family, offset = offset, weights = weights,
                   n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  int.mat <- wrap$int.mat
  
  ybar <- family$linkinv(family$linkfun(y) - offset)
  psi <- c(ybar - muhat + mhat)
  
  if(is.null(span))
    span <- cv_span(a = a, psi = psi, family = family, weights = weights)
  
  out <- sapply(a.vals, loess_est, psi = psi, a = a, span = span, 
                weights = weights, family = family, 
                se.fit = TRUE, int.mat = int.mat)
  
  estimate <- out[1,]
  variance <- out[2,]
  
  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance, span = span)	
  
  return(out)
  
}

# estimate bart outcome model
bart_est <- function(a, y, x, a.vals, weights = NULL, offset = NULL, family = gaussian(),
                     n.iter = 1000, n.adapt = 1000, thin = 10) {
  
  if (is.null(weights))
    weights <- rep(1, nrow(x))
  
  if (is.null(offset))
    offset <- rep(0, nrow(x))
  
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
  
  muhat.mat <- sapply(a, function(a.tmp, ...) {
    
    # for simulations
    xa.tmp <- data.frame(x = x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa)
    return(colMeans(predict(mumod, newdata = xa.tmp, type = "ev")))
    
  })
  
  mhat <- colMeans(muhat.mat)
  
  # integration matrix
  mhat.mat <- matrix(rep(mhat, n), byrow = T, nrow = n)
  int.mat <- muhat.mat - mhat.mat
  
  out <- list(muhat = muhat, mhat = mhat, int.mat = int.mat)
  
  return(out)
  
}

# LOESS function
loess_est <- function(newa, a, psi, family = gaussian(), span, weights = NULL, 
                      se.fit = FALSE, int.mat = NULL, a.vals = NULL, approx = FALSE) {

  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  # subset index
  a.std <- a - newa
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  
  # subset
  a.std <- a.std[idx]
  psi <- psi[idx]
  weights <- weights[idx]
  max.a.std <- max(abs(a.std))
  
  # construct kernel weight
  k.std <- weights*c((1 - abs(a.std/max.a.std)^3)^3)
  gh <- cbind(1, a.std)
  
  # optimize
  bh <- optim(par = c(0,0), fn = opt_fun, k.std = k.std, psi = psi, gh = gh, family = family)
  mu <- family$linkinv(c(bh$par[1]))
  
  if (se.fit & !is.null(int.mat)) {
    
    if (approx) {
    
      kern.mat <- matrix(rep(c((1 - abs((a.vals - newa)/max.a.std)^3)^3), k), byrow = T, nrow = k)
      kern.mat[matrix(rep(abs(a.vals - newa)/max.a.std, k), byrow = T, nrow = k) > 1] <- 0
      g.vals <- matrix(rep(c(a.vals - newa), k), byrow = T, nrow = k)
      
      intfn1.mat <- kern.mat * int.mat[idx,]
      intfn2.mat <- g.vals * kern.mat * int.mat[idx,]
      
      int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
      int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
      
      Dh <- solve(t(gh) %*% diag(k.std) %*% gh)
      V <- crossprod(cbind(k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int1,
                           a.std * k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int2))
      
      sig <- Dh%*%V%*%Dh
      
      return(c(mu = mu, sig = sig[1,1]))
    
    } else {
        
      intfn1.mat <- k.std*t(int.mat[idx,idx])
      intfn2.mat <- a.std*k.std*t(int.mat[idx,idx])
      int1 <- colMeans(intfn1.mat)
      int2 <- colMeans(intfn2.mat)
      
      Dh <- solve(t(gh) %*% diag(k.std) %*% gh)
      V <- crossprod(cbind(k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int1,
                           a.std * k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int2))
      
      sig <- Dh%*%V%*%Dh
      
      return(c(mu = mu, sig = sig[1,1]))

    }
    
  } else
    return(mu)
  
}

# k-fold cross validation to select span
cv_span <- function(a, psi, family = gaussian(), folds = 10, weights = NULL) {
  
  n <- length(a)
  fdx <- sample(x = folds, size = min(n, 1000), replace = TRUE)
  idx <- sample(x = n, size = min(n, 1000), replace = FALSE) # for big data
  a.sub <- a[idx]
  psi.sub <- psi[idx]
  
  cv.mat <- sapply(seq(0.05, 1, by = 0.05), function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a.sub[fdx == k], loess_est, psi = psi.sub[fdx != k], a = a.sub[fdx != k], 
                      weights = weights[fdx != k], span = h, family = family, se.fit = FALSE)
      cv.vec[k] <- mean((psi.sub[fdx == k] - preds)^2, na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })

  cv.err <- colMeans(cv.mat)
  span <- span.seq[which.min(cv.err)]
 
  return(span)
  
}

# optimization used in dr_est
opt_fun <- function(par, k.std, psi, gh, family) {
  
  sum(k.std*(psi - family$linkinv(c(gh %*% par)))^2)
  
}
