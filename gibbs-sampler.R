
gibbs_dr <- function(s, s.id, y.id, t = NULL, x = NULL, w = NULL, 
                     shape = 1e-4, rate = 1e-4, scale = 1e4, 
                     thin = 10, n.iter = 10000, n.adapt = 1000,
                     sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "SL.gam")) {
  
  if (!is.null(w) & !is.null(t)) {
    
    # set up evaluation points & matrices for predictions
    ws <- data.frame(w, s)
    ws.tmp <- ws[!is.na(t),]
    t.tmp <- t[!is.na(t)]
    
    # estimate nuisance outcome model with SuperLearner
    mumod <- SuperLearner(Y = t.tmp, X = ws.tmp, SL.library = sl.lib)
    s <- c(predict(mumod, newdata = ws)$pred)
    s[!is.na(t)] <- t[!is.na(t)]
    
  } else if (!is.null(w) & is.null(t)) {
    
    stop("!is.null(w) & is.null(t)")
    
  } else if (is.null(w) & !is.null(t)) {
    
    s[!is.na(t)] <- t[!is.na(t)]
    warning("replaced values of s with t wherever available.")
    
  }
  
  # remove any s.id not present in y.id
  id <- unique(y.id)[order(unique(y.id))]
  su.id <- unique(s.id)[order(unique(s.id))]
  su.id <- su.id[su.id %in% id]
  
  if(length(su.id) != length(id))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
  
  if(!all(su.id == id))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
    
  s <- s[s.id %in% id]
  s.id <- s.id[s.id %in% id]
  
  if (!is.null(x)) {
    
    d <- aggregate(model.matrix(~ ., data = data.frame(x)), by = list(y.id), mean)
    design <- as.matrix(d[,2:ncol(d)])
    
  } else {
    
    design <- matrix(rep(1, length(id)), ncol = 1)
    
  }
  
  # dimensions
  l <- length(s.id)
  m <- length(id)
  n <- length(y.id)
  p <- ncol(design)

  # initialize exposures
  a <- aggregate(s, by = list(s.id), mean)[,2]
  a_s <- rep(NA, length(s.id))
  
  for (g in id)
    a_s[s.id == g] <- a[id == g]
  
  # initialize parameters
  tau2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2[1] <- var(s)
  beta[1,] <- coef(lm(a ~ 0 + design))
  sigma2[1] <- sigma(lm(a ~ 0 + design))^2
  
  amat <- piHat <- piSig <- matrix(NA, nrow = n.iter + n.adapt, ncol = m)
  amat_y <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  amat_s <- matrix(NA, nrow = n.iter + n.adapt, ncol = l)

  # gibbs sampler for predictors
  for(j in 2:(n.iter + n.adapt)) {
    
    a <- amat[j,] <- sapply(id, function(g, ...) {
      
      sig <- sqrt((sum(s.id == g)/tau2[j - 1] + 1/sigma2[j - 1])^(-1))
      hat <- (sum(design[id == g,]*beta[j - 1,])/sigma2[j - 1] + sum(s[s.id == g])/tau2[j - 1])*sig^2
      return(rnorm(1, hat, sig))
      
    })
    
    for (g in id) {
      
      a_s[s.id == g] <- amat_s[j,s.id == g] <- a[id == g]
      amat_y[j,y.id == g] <- amat[j,id == g] 
      
    }
    
    # Sample EPE params
    tau2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - a_s)^2)/2)
    
    # Sample GPS params
    beta_var <- solve(t(design)%*%design + diag(sigma2[j - 1]/scale, p, p))
    beta[j,] <- rmvnorm(1, beta_var%*%t(design)%*%a, sigma2[j - 1]*beta_var)
    
    sigma2[j] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum(c(a - design%*%beta[j,])^2)/2)

  }
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  amat <- amat[keep,]
  amat_y <- amat_y[keep,]
  amat_s <- amat_s[keep,]
  
  mcmc <- list(beta = beta, sigma2 = sigma2, tau2 = tau2, amat = amat, amat_y = amat_y, amat_s = amat_s)
  
  return(mcmc)
  
}
