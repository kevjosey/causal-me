gp_dr <- function(s, x, s.id, y.id, fmla, shape = 1e-4, rate = 1e-4, scale = 1e4, thin = 10, n.iter = 10000, n.adapt = 1000) {
  
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
  
  d_tmp <- aggregate(model.matrix(as.formula(fmla), data = data.frame(x)), by = list(y.id), mean)
  d <- as.matrix(d_tmp[,2:ncol(d_tmp)])
  
  whichCat <- apply(d, 2, function(x, ...) length(unique(x)) <= 2)
  
  p <- ncol(d)
  pCont <- sum(!whichCat)
  pCat <- sum(whichCat)
  xCont <- d[,!whichCat, drop = FALSE]
  design <- d[,whichCat, drop = FALSE]
  
  if (pCont == 0)
    stop("GP model should not be used with zero continuous covariates")
  
  # dimensions
  l <- length(s.id)
  m <- length(id)
  n <- length(y.id)
  
  # initialize exposures
  a <- aggregate(s, by = list(s.id), mean)[,2]
  a_s <- rep(NA, length(s.id))
  
  for (g in id)
    a_s[s.id == g] <- a[id == g]
  
  # initialize parameters
  f <- array(NA, dim = c(n.iter + n.adapt, m, pCont))
  tau2 <- rep(NA, n.adapt + n.iter)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = pCat)
  sigma2 <- rep(NA, n.adapt + n.iter)
  
  ## starting values
  f[1,,] <- 0
  beta[1,] <- rep(0, pCat)
  sigma2[1] <- 1
  tau2[1] <- 1
  
  for (j in 1 : pCont)
    if (length(unique(xCont[,j])) < m) xCont[,j] <- xCont[,j] + rnorm(m, sd=.00001)
  
  ## We need to do SVD outside of MCMC first
  kernMat <- array(NA, dim = c(pCont,m,m))
  vecMat <- array(NA, dim = c(pCont,m,m))
  invMat <- array(NA, dim = c(pCont,m,m))
  valMat <- array(NA, dim = c(pCont,m))
  
  for (j in 1 : pCont) {
    
    kernMat[j,,] <- diag(m)
    
    for (m1 in 1 : m) {
      for (m2 in 1 : m1) {
        kernMat[j,m1,m2] <- exp(-abs(xCont[m1,j] - xCont[m2,j])/1)
        kernMat[j,m2,m1] <- kernMat[j,m1,m2]
      }
    }
    
    eig <- eigen(kernMat[j,,])
    valMat[j,] <- eig$values
    vecMat[j,,] <- eig$vectors
    invMat[j,,] <- solve(kernMat[j,,])
    
  }
  
  amat <- piHat <- piSig <- matrix(NA, nrow = n.iter + n.adapt, ncol = m)
  amat_y <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  
  for (k in 2 : (n.iter + n.adapt)) {
    
    a <- amat[k,] <- sapply(id, function(g, ...) {
      
      sig <- sqrt((sum(s.id == g)/tau2[k - 1] + 1/sigma2[k - 1])^(-1))
      hat <- ((sum(design[id == g,]*beta[k-1,]) + sum(f[k-1,id == g,]))/sigma2[k - 1] + sum(s[s.id == g])/tau2[k - 1])*sig^2
      return(rnorm(1, hat, sig))
      
    })
    
    for (g in id) {
      
      a_s[s.id == g] <- a[id == g]
      amat_y[k,y.id == g] <- amat[k,id == g] 
      
    }
    
    ## EPE model
    
    tau2[k] <- 1/rgamma(1, shape = shape + l/2, rate = rate + sum((s - a_s)^2)/2)
    
    ## GPS model
    
    # update sigma
    shapeA <- shape + m/2
    rateA <- rate + sum(c(a - design %*% beta[k-1,] - rowSums(f[k-1,,]))^2)/2

    sigma2[k] <- 1/rgamma(1, shapeA, rateA)
    
    ## update categorical params
    aStar <- a - rowSums(f[k - 1,,])
    priorA <- diag(sigma2[k - 1]/scale, pCat, pCat)
    
    beta_var <- solve(t(design)%*%design + priorA)
    beta[k,] <- rmvnorm(1, beta_var%*%t(design)%*%aStar, sigma2[k - 1]*beta_var)
    
    ## update f(x_j)
    f[k,,] <- fUpdate(a = a, fmat = f[k-1,,], design = design, par = beta[k,],
                      sigma2 = sigma2[k], valMat = valMat, vecMat = vecMat)
    
    # f[k,,] <- f[k-1,,]
    # 
    # for (j in 1:ncol(f[k,,])) {
    #   
    #   fj <- as.matrix(f[k,,-j])
    #   rfsum <- rowSums(fj)
    #   star <- a - design %*% beta[k,] - rfsum
    #   temp <- valMat[j,]/(1 + valMat[j,]);
    #   err <- c(rmvnorm(1, rep(0, length(temp)), diag(sqrt(sigma2[k]*temp))))
    #   
    #   # t(t(vecMat[j,,]) * tempObj1) %*% t(vecMat[j,,])
    #   varMat <- vecMat[j,,] %*% diag(temp) %*% t(vecMat[j,,])
    #   f[k,,j] <- varMat %*% star + vecMat[j,,] %*% err
    #   
    # }
    
  }
  
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  f <- f[keep,,]
  amat <- amat[keep,]
  amat_y <- amat_y[keep,]
  
  mcmc <- list(beta = beta, sigma2 = sigma2, tau2 = tau2, f = f, amat = amat, amat_y = amat_y)
  
  return(mcmc)
  
}

