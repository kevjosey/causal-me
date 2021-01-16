
gibbs_dr <- function(s, y, x, s.id, y.id, fmla.s, fmla.a,
                     a.vals = seq(min(a), max(a), length.out = 20), 
                     shape = 1e-4, rate = 1e-4, scale = 1e4, 
                     thin = 10, n.iter = 10000, n.adapt = 1000,
                     n.boot = 100, degree = 2, lambda = seq(0.1, 2.1, by = 0.25),
                     span = NULL, span.seq = seq(0.5, 1, by = 0.05), k = 10,
                     sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "sL.gam")) {
  
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
  
  # initialize exposures
  z <- aggregate(s, by = list(s.id), mean)[,2]
  z_y <- rep(NA, length(y.id))
  z_s <- rep(NA, length(s.id))
  
  for (g in id)
    z_s[s.id == g] <- z[id == g]
  
  # initialize parameters
  tau2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  sigma2 <- rep(NA, n.adapt + n.iter)
  
  beta[1,] <- coef(lm(z ~ 0 + designZ))
  sigma2[1] <- sigma(lm(z ~ 0 + designZ))^2
  tau2[1] <- var(s)
  
  estPost <- matrix(NA, nrow = n.iter + n.adapt, ncol = length(a.vals))
  accept <- rep(NA, n.iter + n.adapt)
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
    tau2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - z_s)^2)/2)
    
    # Sample GPS params
    beta_var <- solve(t(designZ)%*%designZ + diag(sigma2[j - 1]/scale, p, p))
    beta[j,] <- rmvnorm(1, beta_var%*%t(designZ)%*%z, sigma2[j - 1]*beta_var)
    
    sigma2[j] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum(c(z - designZ%*%beta[j,])^2)/2)

  }
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  piSig <- piSig[keep,]
  piHat <- piHat[keep,]
  zMat <- zMat[keep,]
  zMat_y <- zMat_y[keep,]
  
  ### BEFORE PROCEEDING SEE PLUMMER'S SLIDES ON HITTING 
  ### A MOVING TARGET BY Cutting MODEL FEEDBACK
  
  # Using the bayes estimated propensity 
  # estPost <- mclapply.hack(1:nrow(zMat), function(k, ...){
  # 
  #   z <- zMat[k,]
  #   z_y <- zMat_y[k,]
  # 
  #   designY <- data.frame(x, z = z_y)
  #   muMod <- SuperLearner(Y = y, X = designY, SL.library = sl.lib, family = binomial())
  #   muMod_vals <- c(muMod$SL.predict)
  # 
  #   muHat_vals <- sapply(z, function(z_tmp, ...){
  # 
  #     dY <- data.frame(x, z = z_tmp)
  #     colnames(dY) <- colnames(designY)
  #     return(c(predict(muMod, newdata = dY)$pred))
  # 
  #   })
  # 
  #   # aggregate a.vals predictions
  #   dT_tmp <- setDT(data.frame(id = y.id, y = y, muHat = muMod_vals, muHat_mat = muHat_vals))
  #   dT <- dT_tmp[,lapply(.SD, mean), by = y.id][order(id)]
  #   muHat_mat <- dT[,5:ncol(dT)]
  #   mHat_mat <- matrix(rep(colMeans(muHat_mat), m), byrow = TRUE, nrow = m)
  #   int <- rowMeans(muHat_mat - mHat_mat)
  # 
  #   piHat <- c(piHat[k,])
  #   piSig <- c(piSig[k,])
  #   piDens <- dnorm(z, piHat, piSig)
  #   pDens <- sapply(1:length(z), function(i, ...) mean(dnorm(z[i], piHat, piSig)))
  # 
  #   ratio <- (pDens / piDens)
  # 
  #   psi <- ratio*(dT$y - dT$muHat) + mHat_mat[1,]
  #   fit <- sapply(a.vals, dr_est, psi = psi, a = z, int = int, w = size, span = 0.8, se.fit = TRUE)
  # 
  #   return(fit)
  # 
  # })
  
  # estPost <- mclapply.hack(1:nrow(zMat), function(k, ...){
  # 
  #   a_try <- colMeans(zMat_y) + eps_y[k,]
  #   hct_dr(y = y, a = a_try, x = x, y.id = y.id, a.vals = a.vals, span = 0.8, sl.lib = sl.lib)
  # 
  # })
  # 
  # estimate <- t(matrix(unlist(lapply(estPost, function(x) x$estimate)), ncol = length(estPost)))
  # variance <- t(matrix(unlist(lapply(estPost, function(x) x$variance)), ncol = length(estPost)))
  # mcmc <- list(beta = beta, sigma2 = sigma2, tau2 = tau2, estimate = estimate, variance = variance)
  # 
  # out <- list(est = colMeans(estimate),
  #           se1 = apply(estimate, 2, sd),
  #           se2 = sqrt(colMeans(variance)),
  #           ci = apply(estimate, 2, hpd),
  #           mcmc = mcmc)
  
  # estPost <- hct_dr(y = y, a = colMeans(zMat_y), x = x, y.id = y.id, a.vals = a.vals, k = 5, sl.lib = sl.lib)
  
  estPost <- simex_dr(z = colMeans(zMat), y = y, x = x, id = id, y.id = y.id, 
                      sigma = sqrt(mean(sigma2)/table(s.id)), n.boot = n.boot, degree = degree,
                      lambda = lambda, a.vals = a.vals, span = span, span.seq = span.seq, k = k, sl.lib = sl.lib)
  
  estimate <- estPost$estimate
  variance <- estPost$variance

  out <- list(est = estimate, se = sqrt(variance))

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
