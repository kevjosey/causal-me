
gibbs_dr <- function(s, t, w, s.id, y, x, y.id, fmla.s, fmla.z, fmla.y, 
                     a = 1e-4, b = 1e-4, c = 1e4, h = 1, thin = 10,
                     a.vals = seq(quantile(s, .05), quantile(s, .95), length = 20)) {
  
  lvl_y <- unique(y.id)[order(unique(y.id))]
  lvl_s <- unique(s.id)[order(unique(s.id))]
  
  lvl_s <- lvl_s[lvl_s %in% lvl_y]
  lvl <- lvl_y
  
  if(length(lvl_s) != length(lvl))
    stop("some levels in y.id are not represented by s.id. There is no exposure data for these entries.")
  
  if(!all(lvl_s == lvl))
    stop("some levels in y.id are not represented by s.id. There is no exposure data for these entries.")
  
  s <- s[s.id %in% lvl]
  t <- t[s.id %in% lvl]
  w <- w[s.id %in% lvl,]
  s.id <- s.id[s.id %in% lvl]
  
  size <- unname(table(y.id))/mean(table(y.id))
  
  designS <- model.matrix(as.formula(fmla.s), data = data.frame(w))[,-1]
  dZ <- setDT(data.frame(y.id = y.id, model.matrix(as.formula(fmla.z), data = data.frame(x))))
  designZ <- as.matrix(dZ[,lapply(.SD, mean), by = y.id][order(y.id), 2:ncol(dZ)])
  
  l <- nrow(designS)
  m <- nrow(designZ)
  p <- ncol(designS)
  q <- ncol(designZ)

  z <- aggregate(s, by = list(s.id), mean)[,2]

  z_y <- rep(NA, length(y.id))
  z_s <- rep(NA, length(s.id))
  
  for (g in lvl) {
    z_y[y.id == g] <- z[lvl == g]
    z_s[s.id == g] <- z[lvl == g]
  }

  designY <- model.matrix(fmla.y, data = data.frame(x, z = z_y))
  
  n <- nrow(designY)
  r <- ncol(designY)

  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  phi2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  gamma <- matrix(NA, nrow = n.iter + n.adapt, ncol = r)
  
  beta[1,] <- coef(lm(z ~ 0 + designZ))
  sigma2[1] <- sigma(lm(z ~ 0 + designZ))^2
  alpha[1,] <- coef(lm(s ~ 0 + designS))
  phi2[1] <- sigma(lm(s ~ 0 + designS))^2
  gamma[1,] <- coef(glm(y ~ 0 + designY, family = "binomial"))
  
  mu <- plogis(designY%*%gamma[1,])
  Sig <- vcov(glm(y ~ 0 + designY, family = "binomial"))
  
  estPost <- matrix(NA, nrow = n.iter + n.adapt, ncol = length(a.vals))
  zMat <- matrix(NA, nrow = n.iter + n.adapt, ncol = m)
  zMat_y <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)

  # gibbs sampler for predictors
  for(k in 2:(n.iter + n.adapt)) {
  
    # print(k)
    pred <- ifelse(!is.na(t), t, s - designS%*%alpha[k-1,])
    
    z <- zMat[k,] <- sapply(lvl, function(g, ...) {
      
      rnorm(1, (size[lvl == g]*c(designZ[lvl == g,]%*%beta[k-1,])/sigma2[k-1] + 
                  sum(pred[s.id == g])/phi2[k-1]) /
               (sum(s.id == g)/phi2[k-1] + size[lvl == g]/sigma2[k-1]), 
            sqrt((sum(s.id == g)/phi2[k-1] + size[lvl == g]/sigma2[k-1])^(-1)))
      
    })
    
    for (g in lvl) {
      
      z_y[y.id == g] <- zMat_y[k,y.id == g] <- z[lvl == g]
      z_s[s.id == g] <- z[lvl == g]
    
    }
    
    # Sample EPE params
    alpha_var <- solve(t(designS)%*%designS + diag(phi2[k - 1]/c, p, p))
    alpha[k,] <- rmvnorm(1, alpha_var%*%t(designS)%*%(s - z_s), phi2[k - 1]*alpha_var)
    
    phi2[k] <- 1/rgamma(1, shape = a + l/2, rate = b + sum((s - z_s - designS%*%alpha[k,])^2)/2)
    
    # Sample GPS params
    beta_var <- solve(t(designZ)%*%designZ + diag(sigma2[k - 1]/c, q, q))
    beta[k,] <- rmvnorm(1, beta_var%*%t(designZ)%*%z, sigma2[k - 1]*beta_var)
    
    sigma2[k] <- 1/rgamma(1, shape = a + m/2, rate = b + sum((z - designZ%*%beta[k,])^2)/2)
    
    # Sample OM paramaters
    # designY <- model.matrix(fmla.y, data = data.frame(x, z = z_y))
    
    # mh for binary outcome
    # gamma_ <- rmvnorm(1, gamma[1,], h*Sig)
    # mu_ <- plogis(c(designY %*% gamma_))
    # 
    # log.eps <- sum(dbinom(y, 1, mu_, log = TRUE)) -
    #   sum(dbinom(y, 1, mu, log = TRUE)) + 
    #   dmvnorm(gamma_, gamma[1,], h*Sig, log = TRUE) - 
    #   dmvnorm(gamma[k - 1,], gamma[1,], h*Sig, log = TRUE)
    # 
    # log.eps <- log.eps + dmvnorm(gamma_, rep(0, r), diag(1e+4, r), log = TRUE) - 
    #   dmvnorm(gamma[k - 1,], rep(0, r), diag(1e+4, r), log = TRUE)
    # 
    # if ((log(runif(1)) <= log.eps) & !is.na(log.eps)) {
    #   
    #   gamma[k,] <- gamma_
    #   mu <- mu_
    #   
    # } else {
    #   
    #   gamma[k,] <- gamma[k - 1,]
    #   
    # }
    
    # iterated mh
    gamma_ <- gamma[k - 1,]

    for(j in 1:r) {

      gamma_[j] <- rnorm(1, gamma[1,j], sqrt(h*Sig[j,j]))
      mu_ <- plogis(c(designY %*% gamma_))

      log.eps <- sum(dbinom(y, 1, mu_, log = TRUE)) - 
        sum(dbinom(y, 1, mu, log = TRUE)) + 
        dnorm(gamma_[j], gamma[1,j], sqrt(h*Sig[j,j]), log = TRUE) - 
        dnorm(gamma[k - 1,j], gamma[1,j], sqrt(h*Sig[j,j]), log = TRUE)
      
      log.eps <- log.eps + dnorm(gamma_[j], 0, 1e+4, log = TRUE) - 
        dnorm(gamma[k - 1,j], 0, 1e+4, log = TRUE)

      if ((log(runif(1)) <= log.eps) & !is.na(log.eps)) {

        gamma[k,j] <- gamma_[j]
        mu <- mu_

      } else {

        gamma[k,j] <- gamma[k - 1,j]
        gamma_[j] <- gamma[k - 1,j]

      }
    }
    
  }
  
  accept <- apply(gamma, 2, function(x) mean(diff(x) != 0))
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  alpha <- alpha[keep,]
  phi2 <- phi2[keep]
  gamma <- gamma[keep,]
  zMat <- zMat[keep,]
  zMat_y <- zMat_y[keep,]
  
  # formal estimator
  estPost <- t(sapply(1:length(sigma2), function(k, ...){
    
    z <- zMat[k,]
    z_y <- zMat_y[k,]
    
    piHat <- c(designZ%*%beta[k,])
    piSig <- sqrt(sigma2[k])
    piDens <- dnorm(z, piHat, piSig)
    dY <- model.matrix(fmla.y, data = data.frame(x, z = z_y))
    
    pDens <- sapply(z, function(z_tmp, ...) mean(dnorm(z_tmp, piHat, piSig)))
    
    muHat_vals <- sapply(z, function(z_tmp, ...) 
      model.matrix(fmla.y, data = data.frame(x, z = z_tmp))%*%gamma[k,])
    
    # aggregate a.vals predictions
    dT_tmp <- setDT(data.frame(y.id = y.id, y = y, muHat = c(dY%*%gamma[k,]), muHat_mat = muHat_vals))
    dT <- dT_tmp[,lapply(.SD, mean), by = y.id][order(y.id)]
    
    muHat_mat <- dT[,4:ncol(dT)]
    mHat_mat <- matrix(rep(colMeans(muHat_mat), m), byrow = TRUE, nrow = m)
    int <- rowMeans(muHat_mat - mHat_mat)
    
    lambda <- unname(table(y.id))
    wts <- lambda/mean(lambda)
    ratio <- (pDens / piDens)
    
    psi <- ratio*(dT$y - dT$muHat) + mHat_mat[1,]
    fit <- sapply(a.vals, dr_est, psi = psi, a = z, int = int, w = wts, span = 0.9, se.fit = TRUE)
    
    return(fit) 
    
  }, simplify=FALSE))
  
  estimate <- t(matrix(unlist(lapply(estPost, function(x) x[1,])), ncol = length(estPost)))
  variance <- t(matrix(unlist(lapply(estPost, function(x) x[2,])), ncol = length(estPost)))
  
  mcmc <- list(beta = beta, sigma2 = sigma2, alpha = alpha, phi2 = phi2, gamma = gamma, 
               estimate = estimate, variance = variance)
    
  out <- list(est = colMeans(estimate),
            se1 = apply(estimate, 2, sd),
            se2 = sqrt(colMeans(variance)),
            ci = apply(estimate, 2, hpd),
            mcmc = mcmc, accept = accept)
  
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
