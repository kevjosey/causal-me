
gibbs_dr <- function(s, y, x, w = NULL, s.id, y.id, fmla.s = NULL, fmla.a,
                     shape = 1e-4, rate = 1e-4, scale = 1e4, 
                     thin = 10, n.iter = 10000, n.adapt = 1000) {
  
  # remove any s.id not present in y.id
  id <- unique(y.id)[order(unique(y.id))]
  su.id <- unique(s.id)[order(unique(s.id))]
  su.id <- su.id[su.id %in% id]
  
  if(length(su.id) != length(id))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
  
  if(!all(su.id == id))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
  
  if (!is.null(w) & !is.null(fmla.s)){
    
    designS <- model.matrix(as.formula(fmla.s), data = data.frame(w))[,-1]
    designS <- designS[s.id %in% id]
    q <- ncol(designS)
    
  } else{
    
    designS <- NULL
    
  }
    
  s <- s[s.id %in% id]
  s.id <- s.id[s.id %in% id]
  
  dA <- aggregate(model.matrix(as.formula(fmla.a), data = data.frame(x)), by = list(y.id), mean)
  designA <- as.matrix(dA[,2:ncol(dA)])
  
  # dimensions
  l <- length(s)
  m <- length(id)
  n <- length(y)
  p <- ncol(designA)

  # initialize exposures
  a <- aggregate(s, by = list(s.id), mean)[,2]
  a_s <- rep(NA, length(s.id))
  
  for (g in id)
    a_s[s.id == g] <- a[id == g]
  
  # initialize parameters
  tau2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)

  sigma2 <- rep(NA, n.adapt + n.iter)
  
  beta[1,] <- coef(lm(a ~ 0 + designA))
  sigma2[1] <- sigma(lm(a ~ 0 + designA))^2
  
  if (!is.null(designS)) {
   
    alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
    alpha[1,] <- coef(lm(s ~ 0 + designS))
    tau2[1] <- sigma(lm(s ~ 0 + designS))^2
     
  } else {
    
    alpha <- NULL
    tau2[1] <- var(s)
    
  }
  
  amat <- piHat <- piSig <- matrix(NA, nrow = n.iter + n.adapt, ncol = m)
  amat_y <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)

  # gibbs sampler for predictors
  for(j in 2:(n.iter + n.adapt)) {
    
    a <- amat[j,] <- sapply(id, function(g, ...) {
      
      sig <- sqrt((sum(s.id == g)/tau2[j - 1] + 1/sigma2[j - 1])^(-1))
      hat <- (sum(designA[id == g,]*beta[j - 1,])/sigma2[j - 1] + sum(s[s.id == g])/tau2[j - 1])*sig^2
      return(rnorm(1, hat, sig))
      
    })
    
    for (g in id) {
      
      a_s[s.id == g] <- a[id == g]
      amat_y[j,y.id == g] <- amat[j,id == g] 
      
    }
    
    # Sample EPE params
    if (!is.null(designS)) {
      
      alpha_var <- solve(t(designS)%*%designS + diag(tau2[j - 1]/scale, q, q))
      alpha[j,] <- rmvnorm(1, alpha_var%*%t(designS)%*%(s - a_s), tau2[j - 1]*alpha_var)
      
      tau2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - a_s - designS%*%alpha[j,])^2)/2)
    
    } else {
      
      tau2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - a_s)^2)/2)
      
    }
    
    # Sample GPS params
    beta_var <- solve(t(designA)%*%designA + diag(sigma2[j - 1]/scale, p, p))
    beta[j,] <- rmvnorm(1, beta_var%*%t(designA)%*%a, sigma2[j - 1]*beta_var)
    
    sigma2[j] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum(c(a - designA%*%beta[j,])^2)/2)

  }
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  alpha <- alpha[keep,]
  tau2 <- tau2[keep]
  amat <- amat[keep,]
  amat_y <- amat_y[keep,]
  
  mcmc <- list(beta = beta, sigma2 = sigma2, alpha = alpha, tau2 = tau2, amat = amat, amat_y = amat_y)
  
  return(mcmc)
  
}
