
gibbs_dr <- function(s, s.id, id, t = NULL, x = NULL, w = NULL, 
                     shape = 1e-4, rate = 1e-4, scale = 1e4, 
                     thin = 10, n.iter = 10000, n.adapt = 1000,
                     sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "SL.gam")) {
  
  # if (!is.null(w) & !is.null(t)) {
  #   
  #   # set up evaluation points & matrices for predictions
  #   ws <- data.frame(w, s)
  #   ws.tmp <- ws[!is.na(t),]
  #   t.tmp <- t[!is.na(t)]
  #   
  #   # estimate nuisance outcome model with SuperLearner
  #   mumod <- SuperLearner(Y = t.tmp, X = ws.tmp, SL.library = sl.lib)
  #   s <- c(predict(mumod, newdata = ws)$pred)
  #   s[!is.na(t)] <- t[!is.na(t)]
  #   
  # } else if (!is.null(w) & is.null(t)) {
  #   
  #   stop("!is.null(w) & is.null(t)")
  #   
  # } else if (is.null(w) & !is.null(t)) {
  #   
  #   s[!is.na(t)] <- t[!is.na(t)]
  #   warning("replaced values of s with t wherever available.")
  #   
  # }
  
  s[!is.na(t)] <- t[!is.na(t)]
  
  # remove any s.id not present in id
  su.id <- unique(s.id)[order(unique(s.id))]
  su.id <- su.id[su.id %in% id]
  
  if(length(su.id) != length(id))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
  
  if(!all(su.id == id))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
    
  s <- s[s.id %in% id]
  s.id <- s.id[s.id %in% id]
  
  if (!is.null(x)) {
    
    design <- model.matrix(~ ., data = data.frame(x))
    
  } else {
    
    design <- matrix(rep(1, length(id)), ncol = 1)
    
  }
  
  # dimensions
  m <- length(s.id)
  n <- length(id)
  p <- ncol(design)
  q <- ncol(w)

  # initialize exposures
  a <- aggregate(s, by = list(s.id), mean)[,2]
  a_s <- rep(NA, length(s.id))
  
  for (g in id)
    a_s[s.id == g] <- a[id == g]
  
  # initialize parameters
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  tau2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  sigma2 <- rep(NA, n.adapt + n.iter)
  alpha[1,] <- coef(lm((s - a_s) ~ 0 + w))
  tau2[1] <- sigma(lm((s - a_s) ~ 0 + w))^2
  beta[1,] <- coef(lm(a ~ 0 + design))
  sigma2[1] <- sigma(lm(a ~ 0 + design))^2
  
  amat <- piHat <- piSig <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)

  # gibbs sampler for predictors
  for(j in 2:(n.iter + n.adapt)) {
    
    a <- amat[j,] <- sapply(id, function(g, ...) {
      
      sig <- sqrt((sum(s.id == g)/tau2[j - 1] + 1/sigma2[j - 1])^(-1))
      hat <- (sig^2)*(sum(design[id == g,]*beta[j - 1,])/sigma2[j - 1] + 
                sum(is.na(t[s.id == g])*(s[s.id == g] - w[s.id == g,]%*%alpha[j-1,]))/tau2[j - 1] +
                sum(!is.na(t[s.id == g])*s[s.id == g])/tau2[j - 1])
      return(rnorm(1, hat, sig))
      
    })
    
    for (g in id) {
      
      a_s[s.id == g] <- a[id == g]
      
    }
    
    # Sample EPE params
    alpha_var <- solve(t(w[is.na(t),])%*%w[is.na(t),] + diag(tau2[j - 1]/scale, q, q))
    alpha[j,] <- rmvnorm(1, alpha_var%*%t(w[is.na(t),])%*%(s[is.na(t)] - a_s[is.na(t)]), tau2[j - 1]*alpha_var)
    
    tau2[j] <- 1/rgamma(1, shape = shape + m/2, 
                        rate = rate + sum(is.na(t)*(s - a_s - w%*%alpha[j,])^2 +
                                            !is.na(t)*(s - a_s)^2)/2)
    
    # Sample GPS params
    beta_var <- solve(t(design)%*%design + diag(sigma2[j - 1]/scale, p, p))
    beta[j,] <- rmvnorm(1, beta_var%*%t(design)%*%a, sigma2[j - 1]*beta_var)
    
    sigma2[j] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - design%*%beta[j,])^2)/2)

  }
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  alpha <- alpha[keep,]
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  amat <- amat[keep,]
  
  mcmc <- list(beta = beta, sigma2 = sigma2, tau2 = tau2, amat = amat)
  
  return(mcmc)
  
}
