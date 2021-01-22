
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
  
  s <- s[s.id %in% id]
  s.id <- s.id[s.id %in% id]
  
  dZ <- aggregate(model.matrix(as.formula(fmla.a), data = data.frame(x)), by = list(y.id), mean)
  designZ <- as.matrix(dZ[,2:ncol(dZ)])
  
  # dimensions
  l <- length(s)
  m <- length(id)
  n <- length(y)
  p <- ncol(designZ)
  size <- unname(table(y.id))
  
  if (!is.null(w) & !is.null(fmla.s)){
    
    designS <- model.matrix(as.formula(fmla.s), data = data.frame(w))[,-1]
    q <- ncol(designS)
    
  } else
    designS <- NULL
  
  # initialize exposures
  z <- aggregate(s, by = list(s.id), mean)[,2]
  z_s <- rep(NA, length(s.id))
  
  for (g in id)
    z_s[s.id == g] <- z[id == g]
  
  # initialize parameters
  tau2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)

  sigma2 <- rep(NA, n.adapt + n.iter)
  
  beta[1,] <- coef(lm(z ~ 0 + designZ))
  sigma2[1] <- sigma(lm(z ~ 0 + designZ))^2
  
  if (!is.null(designS)) {
   
    alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
    alpha[1,] <- coef(lm(s ~ 0 + designS))
    tau2[1] <- sigma(lm(s ~ 0 + designS))^2
     
  } else {
    
    alpha <- NULL
    tau2[1] <- var(s)
    
  }
  
  zMat <- piHat <- piSig <- matrix(NA, nrow = n.iter + n.adapt, ncol = m)
  zMat_y <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)

  # gibbs sampler for predictors
  for(j in 2:(n.iter + n.adapt)) {
    
    z <- zMat[j,] <- sapply(id, function(g, ...) {
      
      sig <- sqrt((sum(s.id == g)/tau2[j - 1] + 1/sigma2[j - 1])^(-1))
      hat <- (sum(designZ[id == g,]*beta[j - 1,])/sigma2[j - 1] + sum(s[s.id == g])/tau2[j - 1])*sig^2
      return(rnorm(1, hat, sig))
      
    })
    
    for (g in id) {
      
      z_s[s.id == g] <- z[id == g]
      zMat_y[j,y.id == g] <- zMat[j,id == g] 
      
    }
    
    # Sample EPE params
    if (!is.null(designS)) {
      
      alpha_var <- solve(t(designS)%*%designS + diag(tau2[j - 1]/scale, q, q))
      alpha[j,] <- rmvnorm(1, alpha_var%*%t(designS)%*%(s - z_s), tau2[j - 1]*alpha_var)
      
      tau2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - z_s - designS%*%alpha[j,])^2)/2)
    
    } else {
      
      tau2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - z_s)^2)/2)
      
    }
    
    # Sample GPS params
    beta_var <- solve(t(designZ)%*%designZ + diag(sigma2[j - 1]/scale, p, p))
    beta[j,] <- rmvnorm(1, beta_var%*%t(designZ)%*%z, sigma2[j - 1]*beta_var)
    
    sigma2[j] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum(c(z - designZ%*%beta[j,])^2)/2)

  }
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  alpha <- alpha[keep,]
  tau2 <- tau2[keep]
  zMat <- zMat[keep,]
  zMat_y <- zMat_y[keep,]
  
  mcmc <- list(beta = beta, sigma2 = sigma2, alpha = alpha, tau2 = tau2, zMat = zMat, zMat_y = zMat_y)
  
  return(mcmc)
  
}
