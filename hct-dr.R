
# wrapper function to fit a hierarchical, doubly-robust ERC using LOESS regression on a nonparametric model
hct_dr <- function(a, y, x, y.id = NULL,
                   a.vals = seq(min(a), max(a), length.out = 20), 
                   span = NULL, span.seq = seq(0.15, 1, by = 0.05), k = 10,
                   sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.earth")) {	
  
  if (is.null(y.id))
    stop("y.id must be provided.")

  n <- length(y)
  m <- length(unique(y.id))
  
  wrap <- np_est(y = y, a = a, x = x, y.id = y.id, sl.lib = sl.lib)
  
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  pihat <- wrap$pihat
  phat <- wrap$phat
  int <- wrap$int
  a.agg <- wrap$data$a.agg
  y.agg <- wrap$data$y.agg
  wts <- wrap$data$wts
  psi <- (y.agg - muhat)/(pihat/phat) + mhat
  
  if(is.null(span)) {

    folds <- sample(x = k, size = m, replace = TRUE)
    
    cv.mat <- sapply(span.seq, function(h, ...) {
      
      cv.vec <- rep(NA, k)
      
      for(j in 1:k) {
        
        preds <- sapply(a.agg[folds == j], dr_est, psi = psi[folds != j], a = a.agg[folds != j], 
                        int = int[folds != j], wts = wts[folds != j], span = h, se.fit = FALSE)
        cv.vec[j] <- mean(wts[folds == j]*(psi[folds == j] - preds)^2, na.rm = TRUE)
        # some predictions result in `NA` because of the `x` ranges in each fold
        
      }
      
      return(cv.vec)
      
    })
    
    cv.err <- colMeans(cv.mat)
    span <- span.seq[which.min(cv.err)]
    
  }
  
  dr_out <- sapply(a.vals, dr_est, psi = psi, a = a.agg, int = int, wts = wts, span = span, se.fit = TRUE)
  
  estimate <- dr_out[1,]
  variance <- dr_out[2,]

  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance)	
  
  return(out)
  
}

# LOESS function
dr_est <- function(newa, a, psi, int, wts, span, se.fit = FALSE) {
  
  a.std <- a - newa
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  knn <- rep(0, length(a))
  knn[idx] <- 1
  max.a.std <- max(abs(a.std*knn))
  k.std <- wts*((70/81)*(1 - abs(a.std/max.a.std)^3)^3)*knn
  
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
np_est <- function(a, y, x, y.id, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.gam", "SL.earth")){
    
  if (is.null(y.id))
    stop("y.id must be specified.")
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  m <- length(unique(y.id))
  x <- data.frame(x)
  xa <- data.frame(x, a)

  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = y, X = xa, SL.library = sl.lib, family = binomial())
  mumod.vals <- c(mumod$SL.predict)

  # aggregate data and predictions with data.table
  df.tmp <- setDT(data.frame(y.id = y.id, muhat = mumod.vals, y = y, a = a, wts = rep(1, times = n), x = x))
  df <- df.tmp[,lapply(.SD, mean), by = y.id][order(y.id)]
  muhat <- df$muhat
  y.agg <- df$y
  a.agg <- df$a
  x.agg <- as.matrix(df[,6:ncol(df)])
  wts <- df$wts
  
  # estimate nuisance GPS functions via super learner
  pimod <- glm(a.agg ~ x.agg, weights = wts, family = gaussian())
  pimod.vals <- c(pimod$fitted.values)
  pi2mod.vals <- sigma(pimod)^2
  
  # exposure models
  pihat <- dnorm(a.agg, pimod.vals, sqrt(pi2mod.vals))
  phat <- sapply(a.agg, function(a.tmp, ...) weighted.mean(dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals))), w = wts)
  
  # predict marginal outcomes given a.vals (or a.agg)
  muhat.vals <- sapply(a.agg, function(a.tmp, ...) {
  
    xa.tmp <- data.frame(x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa) 
    return(c(predict(mumod, newdata = xa.tmp)$pred))
    
  })
  
  # aggregate muhat.vals and integrate for influence
  dt <- setDT(data.frame(y.id = y.id, muhat = muhat.vals))
  muhat.mat <- dt[,lapply(.SD, mean), by = y.id][order(y.id)][,2:ncol(dt)]
  mhat <- colMeans(muhat.vals)
  mhat.mat <- matrix(rep(mhat, m), byrow = TRUE, nrow = m)
  int <- rowMeans(muhat.mat - mhat.mat)
  
  out <- list(muhat = muhat, mhat = mhat, pihat = pihat, phat = phat, int = int,
              data = list(y.agg = y.agg, a.agg = a.agg, x.agg = x.agg, wts = wts))
  
  return(out)
  
}

# optimization used in dr_est
opt_fun <- function(par, k.std, psi, gh) {
  
  sum(k.std*(psi - plogis(c(gh %*% par)))^2)
  
}