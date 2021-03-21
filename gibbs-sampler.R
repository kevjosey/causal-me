
gibbs_dr <- function(s, star, s.id, id, w = NULL, x = NULL,
                     shape = 1e-4, rate = 1e-4, scale = 1e4, 
                     thin = 10, n.iter = 10000, n.adapt = 1000) {
  
  # remove any s.id not present in id
  su.id <- unique(s.id)[order(unique(s.id))]
  su.id <- su.id[su.id %in% id]
  
  if(length(su.id) != length(id))
    stop("some observations in y.id are not represented by measurements of s.id.\n
         There are no exposure data for these entries.")
  
  if(!all(su.id == id[order(id)]))
    stop("some observations in y.id are not represented by measurements of s.id.\n
         There are no exposure data for these entries.")
  
  s <- s[s.id %in% id]
  star <- star[s.id %in% id]
  w <- w[s.id %in% id,]
  s.id <- s.id[s.id %in% id]
  
  if (is.null(x)) {
    
    x <- matrix(1, nrow = length(id), ncol = 1)
    
  } else {
    
    x <- model.matrix(~ ., data.frame(x))
    
  }
  
  if (is.null(w)) {
    
    w <- cbind(1, star)
    
  } else {
    
    w <- model.matrix(~ ., data.frame(w))
    
  }
  
  w.tmp <- w[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  
  # dimensions
  m <- length(s.id)
  n <- length(id)
  p <- ncol(x)
  q <- ncol(w)
  l <- sum(!is.na(s))
  
  # initialize exposures
  a <- aggregate(star, by = list(s.id), mean)[,2]
  a_s <- rep(NA, length(s.id))
  
  for (g in id)
    a_s[s.id == g] <- a[id == g]
  
  # initialize parameters
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  sigma2 <- rep(NA, n.adapt + n.iter)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + w.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- var(star - a_s)
  omega2[1] <- sigma(lm(s.tmp ~ 0 + w.tmp))^2
  
  amat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  smat <- matrix(NA, nrow = n.iter + n.adapt, ncol = m)
  
  # gibbs sampler for predictors
  for(j in 2:(n.iter + n.adapt)) {
    
    sig <- sqrt((1/tau2[j - 1] + 1/omega2[j - 1])^(-1))
    hat <- (sig^2)*(a_s/tau2[j - 1] + w%*%alpha[j - 1,]/omega2[j - 1])
    t <- smat[j,] <- rnorm(m, hat, sig)
    t[!is.na(s)] <- s[!is.na(s)]
    
    # set up evaluation points & matrices for predictions

    a <- amat[j,] <- sapply(id, function(g, ...) {
      
      idx <- s.id == g
      sig <- sqrt((sum(idx)/tau2[j - 1] + 1/sigma2[j - 1])^(-1))
      hat <- (sig^2)*(sum(t[idx])/tau2[j - 1] + 
        sum(x[id == g,]*beta[j - 1,])/sigma2[j - 1])
      return(rnorm(1, hat, sig))
      
    })
    
    for (g in id)
      a_s[s.id == g] <- a[id == g]
    
    # Sample pred params
    
    alpha_var <- solve(t(w.tmp)%*%w.tmp + diag(omega2[j - 1]/scale, q, q))
    alpha[j,] <- rmvnorm(1, alpha_var%*%t(w.tmp)%*%s.tmp, omega2[j - 1]*alpha_var)
    
    omega2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                            sum((s.tmp - (w.tmp)%*%alpha[j,])^2)/2)
    
    # Sample agg params
    
    tau2[j] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((t - a_s)^2)/2)
    
    # Sample GPS params
    
    beta_var <- solve(t(x)%*%x + diag(sigma2[j - 1]/scale, p, p))
    beta[j,] <- rmvnorm(1, beta_var%*%t(x)%*%a, sigma2[j - 1]*beta_var)
    
    sigma2[j] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - x%*%beta[j,])^2)/2)
    
  }
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  amat <- amat[keep,]
  smat <- smat[keep,]
  
  mcmc <- list(beta = beta, alpha = alpha, sigma2 = sigma2, tau2 = tau2,
               omega2 = omega2, amat = amat, smat = smat)
  
  return(mcmc)
  
}
