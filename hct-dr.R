
## hcTMLE: function to run full hierarchical TMLE under nonparametric model with continuous treatment

hct_dr <- function(a, y, x, y.id = NULL, a.vals = seq(min(a), max(a), length.out = 20), 
                   span = NULL, span.seq = seq(0.15, 1, by = 0.05), k = 10,
                   sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "sL.gam")) {	
  
  n <- nrow(x)
  if (is.null(y.id)) y.id <- 1:n
  m <- length(unique(y.id))
  
  size <- unname(table(y.id))
  attributes(size) <- NULL
  w <- size/mean(size)
  
  wrap <- np_est(y = y, a = a, x = x, y.id = y.id, sl.lib = sl.lib)
  muhat <- wrap$muhat
  mhat <- wrap$mhat
  pihat <- wrap$pihat
  phat <- wrap$phat
  int <- wrap$int
  a.agg <- wrap$data$a.agg
  y.agg <- wrap$data$y.agg
  psi <- (y.agg - muhat)/(pihat/phat) + mhat
  
  if(is.null(span)) {

    folds <- sample(x = k, size = m, replace = TRUE)
    
    cv.mat <- sapply(span.seq, function(h, ...) {
      
      cv.vec <- rep(NA, k)
      
      for(j in 1:k) {
        
        preds <- sapply(a.agg[folds == j], dr_est, psi = psi[folds != j], a = a.agg[folds != j], 
                        int = int[folds != j], w = w[folds != j], span = h, se.fit = FALSE)
        cv.vec[j] <- mean(w[folds == j]*(psi[folds == j] - preds)^2, na.rm = TRUE)
        # some predictions result in `NA` because of the `x` ranges in each fold
        
      }
      
      return(cv.vec)
      
    })
    
    cv.err <- colMeans(cv.mat)
    span <- span.seq[which.min(cv.err)]
    
  }
  
  dr_out <- sapply(a.vals, dr_est, psi = psi, a = a.agg, int = int, w = w, span = span, se.fit = TRUE)
  
  estimate <- dr_out[1,]
  variance <- dr_out[2,]

  names(estimate) <- names(variance) <- a.vals
  out <- list(estimate = estimate, variance = variance)	
  
  return(out)
  
}

dr_est <- function(newa, a, psi, int, w, span, se.fit = FALSE) {
  
  a.std <- a - newa
  k <- floor(min(span, 1)*length(a))
  idx <- order(abs(a.std))[1:k]
  knn <- rep(0, length(a))
  knn[idx] <- 1
  max.a.std <- max(abs(a.std*knn))
  k.std <- w*((70/81)*(1 - abs(a.std/max.a.std)^3)^3)*knn
  
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

np_est <- function(a, y, x, y.id, sl.lib){
  
  if (is.null(y.id))
    stop("y.id must be specified.")
  
  # set up evaluation points & matrices for predictions
  n <- nrow(x)
  m <- length(unique(y.id))
  x <- data.frame(x)
  xa <- data.frame(x, a)
  
  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = y, X = xa, id = y.id, SL.library = sl.lib, family = binomial())
  mumod.vals <- c(mumod$SL.predict)

  # aggregate data and predictions with data.table
  df.tmp <- setDT(data.frame(y.id = y.id, muhat = mumod.vals, y = y, a = a, x = x))
  df <- df.tmp[,lapply(.SD, mean), by = y.id][order(y.id)]
  muhat <- df$muhat
  y.agg <- df$y
  a.agg <- df$a
  x.agg <- as.data.frame(df[,5:ncol(df)])
  
  # estimate nuisance GPS functions via super learner
  pimod <- SuperLearner(Y = a.agg, X = x.agg, SL.library = sl.lib)
  pimod.vals <- c(pimod$SL.predict)
  pi2mod <- SuperLearner(Y = (a.agg - pimod.vals)^2, X = x.agg, SL.library = sl.lib)
  pi2mod.vals <- c(pi2mod$SL.predict)
  
  if (any(pi2mod.vals < 0)) {
    
    pi2mod <- SuperLearner(Y = (a.agg - pimod.vals)^2, X = x.agg, SL.library = "SL.mean")
    pi2mod.vals <- pi2mod$SL.predict
    
  }
  
  # exposure models
  pihat <- dnorm(a.agg, pimod.vals, sqrt(pi2mod.vals))
  phat <- sapply(a.agg, function(a.tmp, ...) mean(dnorm(a.tmp, pimod.vals, sqrt(pi2mod.vals))))
  
  # predict marginal outcomes given a.vals (or a.agg)
  muhat.vals <- sapply(a.agg, function(a.tmp, ...) {
  
    xa.tmp <- data.frame(x, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa)
    return(c(predict(mumod, newdata = xa.tmp)$pred))
    
  })
  
  # aggregate muhat.vals and integrate for influence
  dt <- setDT(data.frame(y.id = y.id, muhat = muhat.vals))
  muhat.mat <- dt[,lapply(.SD, mean), by = y.id][order(y.id)][,2:ncol(dt)]
  mhat <- colMeans(muhat.mat)
  mhat.mat <- matrix(rep(mhat, m), byrow = TRUE, nrow = m)
  int <- rowMeans(muhat.mat - mhat.mat)
  
  out <- list(muhat = muhat,  mhat = mhat, pihat = pihat, phat = phat, int = int,
              data = list(y.agg = y.agg, a.agg = a.agg, x.agg = x.agg))
  
  return(out)
  
}

opt_fun <- function(par, k.std, psi, gh) {
  
  sum(k.std*(psi - plogis(c(gh %*% par)))^2)
  
}
