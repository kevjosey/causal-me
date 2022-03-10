# wrapper function to fit an ERF with measurement error using LOESS regression on a nonparametric models
erf <- function(a, y, x, family = gaussian(), offset = NULL, weights = NULL,
                a.vals = seq(min(a), max(a), length.out = 100),
                n.iter = 10000, n.adapt = 1000, thin = 10,
                span = NULL, span.seq = seq(0.05, 1, by = 0.05), k = 5){	
  
  
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
  psi <- c((ybar - muhat) + mhat)
  
  if(is.null(span)) {

    idx <- sample(x = n, size = min(n, 1000), replace = FALSE)

    a.sub <- a[idx]
    psi.sub <- psi[idx]

    folds <- sample(x = k, size = min(n, 1000), replace = TRUE)

    cv.mat <- sapply(span.seq, function(h, ...) {

      cv.vec <- rep(NA, k)

      for(j in 1:k) {

        preds <- sapply(j, a.sub, loess_est, psi = psi.sub[folds != j], a = a.sub[folds != j], 
                        span = h, family = gaussian(), se.fit = FALSE)
        cv.vec[j] <- mean((psi.sub[folds == j] - preds)^2, na.rm = TRUE)

      }

      return(cv.vec)

    })

    cv.err <- colMeans(cv.mat)
    span <- span.seq[which.min(cv.err)]

  }
  
  out <- sapply(a.vals, loess_est, psi = psi, a = a, span = span, 
                   family = gaussian(), se.fit = TRUE, int.mat = int.mat)
  
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
    offset <- rep(1, nrow(x))
  
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
  
  mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
  
  # exposure model for integtion
  a.std <- c(c(a, a.vals) - mean(a)) / sd(a)
  dens <- density(a.std[1:n])
  phat.vals <- approx(x = dens$x, y = dens$y, xout = a.std[-(1:n)])$y / sd(a)
  
  # integration matrix
  phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)
  mhat.mat <- matrix(rep(colMeans(muhat.mat), n), byrow = T, nrow = n)
  int.mat <- (muhat.mat - mhat.mat) * phat.mat
  
  out <- list(muhat = muhat, mhat = mhat, int.mat = int.mat)
  
  return(out)
  
}

# LOESS function
loess_est <- function(newa, a, psi, span, family = gaussian(), se.fit = FALSE, int.mat = NULL) {

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

  if (se.fit & !is.null(int.mat)) {
    
    kern.mat <- matrix(rep(c((1 - abs((a.vals - newa)/max.a.std)^3)^3), k), byrow = T, nrow = k)
    kern.mat[matrix(rep(abs(a.vals - newa)/max.a.std, k), byrow = T, nrow = k) > 1] <- 0
    g2 <- matrix(rep(c(a.vals - newa), k), byrow = T, nrow = k)
    intfn1.mat <- kern.mat * int.mat[idx,]
    intfn2.mat <- g2 * kern.mat * int.mat[idx,]
    int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                    (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
    int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
                    (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
    
    Dh <- solve(t(gh) %*% diag(k.std) %*% gh)
    V <- crossprod(cbind(k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int1,
                         a.std * k.std * c(psi - family$linkinv(c(gh%*%bh$par))) + int2))
    
    sig <- Dh%*%V%*%Dh
    
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)

}

# optimization used in dr_est
opt_fun <- function(par, k.std, psi, gh, family) {
  
  sum(k.std*(psi - family$linkinv(c(gh %*% par)))^2)
  
}
