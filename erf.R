# wrapper function to fit a nonparametric ERF with measurement error using 
# kernel-weighted least-squares regression
erf <- function(a, y, x, offset = NULL, bart = TRUE,
                a.vals = seq(min(a), max(a), length.out = 100),
                n.iter = 10000, n.adapt = 1000, thin = 10,
                bw = NULL, bw.seq = seq(0.1, 2, by = 0.1), folds = 5) {	
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  weights <- exp(offset)
  ybar <- exp(log(y) - offset)
  n <- length(a)
  
  if (bart)
    wrap <- bart_est(y = ybar, a = a, x = x, a.vals = a.vals, weights = weights,
                     n.iter = n.iter, n.adapt = n.adapt, thin = thin)
  else
    wrap <- glm_est(y = ybar, a = a, x = x, a.vals = a.vals, weights = weights)
  
  psi <- wrap$psi
  int.mat <- wrap$int.mat
  
  # select bw if null
  if (is.null(bw))
    bw <- cv_bw(a = a, psi = psi, folds = folds, bw.seq = bw.seq)

  # asymptotics
  out <- sapply(a.vals, kern_est, psi = psi, a = a, weights = weights,
                bw = bw, se.fit = TRUE, int.mat = int.mat, a.vals = a.vals)
  
  estimate <- out[1,]
  variance <- out[2,]
  
  # bootstrap
  # boot.idx <- cbind(1:n, replicate(200, sample(x = n, size = n, replace = TRUE)))
  # 
  # out <- apply(boot.idx, 2, function(idx, ...) {
  #   
  #   sapply(a.vals, kern_est, psi = psi[idx], a = a[idx], 
  #          weights = weights[idx], bw = bw, se.fit = FALSE)
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
                     n.iter = 1000, n.adapt = 1000, thin = 10, ...) {
  
  if (is.null(weights))
    weights <- rep(1, nrow(x))
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- data.frame(x, a = a)
  colnames(xa) <- c(colnames(x), "a")
  
  # outcome model
  mumod <- dbarts::bart(y.train = y, x.train = xa, weights = weights, keeptrees = TRUE,
                        ndpost = n.iter, nskip = n.adapt, keepevery = thin, verbose = FALSE)
  muhat <- mumod$yhat.train.mean
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    xa.tmp <- data.frame(x = x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa)
    return(colMeans(predict(mumod, newdata = xa.tmp, type = "ev")))
    
  })
  
  mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = weights)
  mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), x = a)$y
  mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
  
  # pseudo outcome
  psi <- c(y - muhat + mhat)
  
  # integration matrix
  a.std <- (c(a, a.vals) - mean(a))/sd(a)
  dens <- density(a.std[1:n])
  phat.vals <- approx(x = dens$x, y = dens$y, xout = a.std[-(1:n)])$y / sd(a)
  phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)  
  int.mat <- (muhat.mat - mhat.mat)*phat.mat
  
  out <- list(psi = psi, int.mat = int.mat)
  
  return(out)
  
}

# estimate glm outcome model
glm_est <- function(a, y, x, a.vals, weights = NULL, ...) {
  
  if (is.null(weights))
    weights <- rep(1, nrow(x))
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  x <- data.frame(x)
  xa <- cbind(as.matrix(x), a - 10, cos(pi*(a - 6)/4), (a - 10)*x[,1])
  
  # outcome model
  mumod <- glm(y ~ xa, weights = weights, family = quasipoisson())
  muhat <- mumod$fitted.values
  
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    
    xa.tmp <- cbind(1, as.matrix(x), a.tmp - 10, cos(pi*(a.tmp - 6)/4), (a.tmp - 10)*x[,1])
    return(c(exp(xa.tmp %*% mumod$coefficients)))
    
  })
  
  mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = weights)
  mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), x = a)$y
  mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
  
  # pseudo outcome
  psi <- c(y - muhat + mhat)
  
  # integration matrix
  a.std <- (c(a, a.vals) - mean(a))/sd(a)
  dens <- density(a.std[1:n])
  phat.vals <- approx(x = dens$x, y = dens$y, xout = a.std[-(1:n)])$y / sd(a)
  phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)  
  int.mat <- (muhat.mat - mhat.mat)*phat.mat
  
  out <- list(psi = psi, int.mat = int.mat)
  
  return(out)
  
}

# LOESS function
kern_est <- function(a.new, a, psi, bw, weights = NULL, se.fit = FALSE, int.mat = NULL, a.vals = NULL) {
  
  n <- length(a)
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  ## LOESS Kernel
  
  # index
  # a.std <- a - a.new
  # k <- floor(min(bw, 1)*length(a))
  # idx <- order(abs(a.std))[1:k]
  
  # subset
  # a.std <- a.std[idx]
  # psi <- psi[idx]
  # weights <- weights[idx]
  # max.a.std <- max(abs(a.std))
  
  # construct kernel weight
  # k.std <- c((1 - abs(a.std/max.a.std)^3)^3)
  # g.std <- cbind(1, a.std)
  
  ## Gaussian kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  # kernel-weighted least squares
  b <- lm(psi ~ -1 + g.std, weights = k.std*weights)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat)) {
  
    eta <- c(g.std %*% b)
    
    ## LOESS
    # kern.mat <- matrix(rep(c((1 - abs((a.vals - a.new)/max.a.std)^3)^3), k), byrow = T, nrow = k)
    # kern.mat[matrix(rep(abs(a.vals - a.new)/max.a.std, k), byrow = T, nrow = k) > 1] <- 0
    # g.vals <- matrix(rep(c(a.vals - a.new), k), byrow = T, nrow = k)
    # intfn1.mat <- kern.mat * int.mat[idx,]
    # intfn2.mat <- g.vals * kern.mat * int.mat[idx,]
    # 
    # int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
    #                 (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
    # int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
    #                 (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
    
    ## Gaussian kernel
    kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
    g.vals <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.vals * kern.mat * int.mat

    int1 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2)
    int2 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                    (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights * (k.std * (psi - eta) + int1),
               weights * (a.std * k.std * (psi - eta) + int2))
    sig2 <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = sig2[1,1]))
    
  } else
    return(mu)
  
}

# k-fold cross validation to select bw
cv_bw <- function(a, psi, weights = NULL, folds = 5, bw.seq = seq(0.05, 1, by = 0.05)) {
  
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
  n <- length(a)
  fdx <- sample(x = folds, size = min(n, 1000), replace = TRUE)
  idx <- sample(x = n, size = min(n, 1000), replace = FALSE) # for big data
  a.sub <- a[idx]
  psi.sub <- psi[idx]
  
  cv.mat <- sapply(bw.seq, function(h, ...) {
    
    cv.vec <- rep(NA, folds)
    
    for(k in 1:folds) {
      
      preds <- sapply(a.sub[fdx == k], kern_est, psi = psi.sub[fdx != k], a = a.sub[fdx != k], 
                      weights = weights[fdx != k], bw = h, se.fit = FALSE)
      cv.vec[k] <- mean((psi.sub[fdx == k] - preds)^2, na.rm = TRUE)
      
    }
    
    return(cv.vec)
    
  })
  
  cv.err <- colMeans(cv.mat)
  bw <- bw.seq[which.min(cv.err)]
  
  return(bw)
  
}

