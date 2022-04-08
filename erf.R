# wrapper function to fit an ERF with measurement error using LOESS regression on a nonparametric models
erf <- function(a, y, x, offset = NULL,
                a.vals = seq(min(a), max(a), length.out = 100),
                n.iter = 10000, n.adapt = 1000, thin = 10,
                span = NULL, span.seq = seq(0.1, 2, by = 0.1), folds = 5) {	
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  weights <- exp(offset)
  ybar <- exp(log(y) - offset)
  n <- length(a)
  
  wrap <- bart_est(y = ybar, a = a, x = x, a.vals = a.vals, weights = weights,
                   n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  
  psi <- wrap$psi
  int.mat <- wrap$int.mat
  
  # select span if null
  if (is.null(span))
    span <- cv_span(a = a, psi = psi, folds = folds, span.seq = span.seq)

  # asymptotics
  out <- sapply(a.vals, loess_est, psi = psi, a = a, weights = weights,
                span = span, se.fit = TRUE, int.mat = int.mat)
  
  estimate <- out[1,]
  variance <- out[2,]
  
  # bootstrap
  # boot.idx <- cbind(1:n, replicate(200, sample(x = n, size = n, replace = TRUE)))
  # 
  # out <- apply(boot.idx, 2, function(idx, ...) {
  #   
  #   sapply(a.vals, loess_est, psi = psi[idx], a = a[idx], 
  #          weights = weights[idx], span = span, se.fit = FALSE)
  #   
  # })
  # 
  # estimate <- out[,1]
  # variance <- apply(out[,2:ncol(out)], 1, var)
  
  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance)	
  
  return(out)
  
}

# estimate bart outcome model
bart_est <- function(a, y, x, a.vals, weights = NULL,
                     n.iter = 1000, n.adapt = 1000, thin = 10) {
  
  if (is.null(weights))
    weights <- rep(0, nrow(x))
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a = a)
  colnames(xa) <- c(colnames(x), "a")
  
  # outcome model
  mumod <- dbarts::bart(y.train = y, x.train = xa, weights = weights, keeptrees = TRUE,
                        ndpost = n.iter, nskip = n.adapt, keepevery = thin, verbose = FALSE)
  muhat <- mumod$yhat.train.mean
  
  muhat.mat <- sapply(a, function(a.tmp, ...) {
    
    xa.tmp <- data.frame(x = x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa)
    return(colMeans(predict(mumod, newdata = xa.tmp, type = "ev")))
    
  })
  
  mhat <- apply(muhat.mat, 2, weighted.mean, w = weights)
  
  # pseudo outcome
  psi <- c(y - muhat + mhat)
  
  # integration matrix
  int.mat <- muhat.mat - matrix(rep(mhat, n), byrow = T, nrow = n)
  
  out <- list(psi = psi, int.mat = int.mat)
  
  return(out)
  
}

# LOESS function
loess_est <- function(a.new, a, psi, span, weights = NULL, se.fit = FALSE, int.mat = NULL) {
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  # standardize
  a.std <- a - a.new
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  
  # subset
  a.std <- a.std[idx]
  psi <- psi[idx]
  weights <- weights[idx]
  
  # construct kernel
  max.a.std <- max(abs(a.std))
  k.std <- c((1 - abs(a.std/max.a.std)^3)^3)
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ -1 + g.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat)) {
    
    eta <- c(g.std %*% b)
    
    int1 <- colMeans(k.std * t(int.mat[idx,idx]))
    int2 <- colMeans(a.std * k.std * t(int.mat[idx,idx]))
    
    U <- solve(t(g.std) %*% diag(weights*k.std) %*% g.std) 
    V <- cbind(weights * (k.std * (psi - eta) + int1),
               weights * (a.std * k.std * (psi - eta) + int2))
    sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)
  
}

# k-fold cross validation to select span
cv_span <- function(a, psi, weights = NULL, folds = 5, span.seq = seq(0.05, 1, by = 0.05)) {
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  n <- length(a)
  fdx <- sample(x = folds, size = min(n, 1000), replace = TRUE)
  idx <- sample(x = n, size = min(n, 1000), replace = FALSE) # for big data
  a.sub <- a[idx]
  psi.sub <- psi[idx]
  
  cv.mat <- sapply(span.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a.sub[fdx == k], loess_est, psi = psi.sub[fdx != k], a = a.sub[fdx != k], 
                      weights = weights[fdx != k], span = h, se.fit = FALSE)
      cv.vec[k] <- mean((psi.sub[fdx == k] - preds)^2, na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  span <- span.seq[which.min(cv.err)]
  
  return(span)
  
}

