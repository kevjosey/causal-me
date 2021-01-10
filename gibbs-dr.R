
gibbs_dr <- function(s, y, x, s.id, y.id, fmla.s, fmla.a, 
                     shape = 1e-4, rate = 1e-4, scale = 1e4, thin = 10, n.iter = 10000, n.adapt = 1000,
                     a.vals = seq(quantile(s, .05), quantile(s, .95), length = 20)) {
  
  # remove any s.id not present in y.id
  lvl <- unique(y.id)[order(unique(y.id))]
  slvl <- unique(s.id)[order(unique(s.id))]
  slvl <- slvl[slvl %in% lvl]
  
  if(length(slvl) != length(lvl))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
  
  if(!all(slvl == lvl))
    stop("some observations in y.id are not represented by measurements of s.id. There is no exposure data for these entries.")
  
  s <- s[s.id %in% lvl]
  s.id <- s.id[s.id %in% lvl]
  size <- unname(table(y.id))
  
  dZ <- aggregate(model.matrix(as.formula(fmla.a), data = data.frame(x)), by = list(y.id), mean)
  designZ <- as.matrix(dZ[,2:ncol(dZ)])
                  
  l <- length(s)
  m <- length(size)
  n <- length(y)
  p <- ncol(designZ)

  z <- aggregate(s, by = list(s.id), mean)[,2]
  z_y <- rep(NA, length(y.id))
  z_s <- rep(NA, length(s.id))
  
  for (g in lvl) {
    z_y[y.id == g] <- z[lvl == g]
    z_s[s.id == g] <- z[lvl == g]
  }
  
  tau2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  sigma2 <- rep(NA, n.adapt + n.iter)
  
  beta[1,] <- coef(lm(z ~ 0 + designZ))
  sigma2[1] <- sigma(lm(z ~ 0 + designZ))^2
  tau2[1] <- var(s)
  
  estPost <- matrix(NA, nrow = n.iter + n.adapt, ncol = length(a.vals))
  accept <- rep(NA, n.iter + n.adapt)
  zMat <- matrix(NA, nrow = n.iter + n.adapt, ncol = m)
  zMat_y <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)

  # gibbs sampler for predictors
  for(k in 2:(n.iter + n.adapt)) {
    
    z <- zMat[k,] <- sapply(lvl, function(g, ...) {
      
      rnorm(1, (c(designZ%*%beta[k-1,])/sigma2[k-1] + sum(s[s.id == g])/tau2[k-1]) / 
              (sum(s.id == g)/tau2[k-1] + 1/sigma2[k-1]), 
            sqrt((sum(s.id == g)/tau2[k-1] + 1/sigma2[k-1])^(-1)))
      
    })
    
    for (g in lvl) {
      
      z_y[y.id == g] <- zMat_y[k,y.id == g] <- z[lvl == g]
      z_s[s.id == g] <- z[lvl == g]
    
    }
    
    # Sample EPE params
    tau2[k] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - z_s)^2)/2)
    
    # Sample GPS params
    beta_var <- solve(t(designZ)%*%designZ + diag(sigma2[k - 1]/scale, p, p))
    beta[k,] <- rmvnorm(1, beta_var%*%t(designZ)%*%z, sigma2[k - 1]*beta_var)
    
    sigma2[k] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum(c(z - designZ%*%beta[k,])^2)/2)
    
  }
  
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  zMat <- zMat[keep,]
  zMat_y <- zMat_y[keep,]
  
  SL.library <- c("SL.mean", "SL.glm", "SL.glm.interaction")
  
  ### BEFORE PROCEEDING SEE PLUMMER'S SLIDES ON HITTING 
  ### A MOVING TARGET BY Cutting MODEL FEEDBACK
  
  # Using the bayes estimated propensity 
  estPost <- mclapply.hack(1:nrow(zMat), function(k, ...){

    z <- zMat[k,]
    z_y <- zMat_y[k,]
    
    designY <- data.frame(x, z = z_y)
    muMod <- SuperLearner(Y = y, X = designY, id = y.id, SL.library = SL.library, family = binomial())
    muMod_vals <- c(muMod$SL.predict)
    
    # Sample OM paramaters
    # for (j in 2:(n.iter + n.adapt)) {
    #   
    #   gamma_ <- c(rmvnorm(1, c(gamma_tmp[j - 1,]), h*Sig))
    #   mu_ <- plogis(c(designY %*% gamma_))
    #   
    #   log.eps <- sum(dbinom(y, 1, mu_, log = TRUE)) -
    #     sum(dbinom(y, 1, mu, log = TRUE))
    #   
    #   log.eps <- log.eps + dmvnorm(gamma_, rep(0, q), diag(1e+4, q), log = TRUE) -
    #     dmvnorm(gamma_tmp[j - 1,], rep(0, q), diag(1e+4, q), log = TRUE)
    #   
    #   if ((log(runif(1)) <= log.eps) & !is.na(log.eps)) {
    #     
    #     gamma_tmp[j,] <- gamma_
    #     mu <- mu_
    #     
    #   } else {
    #     
    #     gamma_tmp[j,] <- gamma_tmp[j - 1,]
    #     
    #   }
    #   
    # }
    # 
    # accept <- mean(diff(gamma_tmp[,1]) != 0)
    # gamma <- colMeans(gamma_tmp[keep,])
    # muMod_vals <- plogis(c(designY%*%gamma))
    
    muHat_vals <- sapply(z, function(z_tmp, ...){
      
      dY <- data.frame(x, z = z_tmp)
      colnames(dY) <- colnames(designY)
      return(c(predict(muMod, newdata = dY)$pred))
      
    })
    
    # aggregate a.vals predictions
    dT_tmp <- setDT(data.frame(id = y.id, y = y, muHat = muMod_vals, muHat_mat = muHat_vals))
    dT <- dT_tmp[,lapply(.SD, mean), by = y.id][order(id)]
    muHat_mat <- dT[,5:ncol(dT)]
    mHat_mat <- matrix(rep(colMeans(muHat_mat), m), byrow = TRUE, nrow = m)
    int <- rowMeans(muHat_mat - mHat_mat)
    
    piHat <- c(designZ %*% beta[k,])
    piSig <- c(sigma2[k])
    piDens <- dnorm(z, piHat, piSig)
    pDens <- sapply(1:length(z), function(i, ...) mean(dnorm(z[i], piHat, piSig)))

    ratio <- (pDens / piDens)

    psi <- ratio*(dT$y - dT$muHat) + mHat_mat[1,]
    fit <- sapply(a.vals, dr_est, psi = psi, a = z, int = int, w = size, span = 0.8, se.fit = TRUE)

    return(fit)

  })
  
  estimate <- t(matrix(unlist(lapply(estPost, function(x) x[1,])), ncol = length(estPost)))
  variance <- t(matrix(unlist(lapply(estPost, function(x) x[2,])), ncol = length(estPost)))
  # accept <- c(unlist(lapply(estPost, function(x) x$accept)), ncol = length(estPost))
  
  mcmc <- list(beta = beta, sigma2 = sigma2, tau2 = tau2, 
               estimate = estimate, variance = variance)
    
  out <- list(est = colMeans(estimate),
            se1 = apply(estimate, 2, sd),
            se2 = sqrt(colMeans(variance)),
            ci = apply(estimate, 2, hpd),
            mcmc = mcmc)
  
  return(out)
  
}

# Highest Posterior Density
hpd <- function(x, alpha = 0.05){
  
  n <- length(x)
  m <- round(n * alpha)
  x <- sort(x)
  y <- x[(n - m + 1):n] - x[1:m]
  z <- min(y)
  k <- which(y == z)[1]
  c(x[k], x[n - m + k])
  
}
