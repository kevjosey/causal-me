
blp <- function(s, s.id, x = NULL) {
  
  wts <- c(unname(table(s.id)))
  
  # initialize exposures
  z_tmp <- aggregate(s, by = list(s.id), mean)
  id <- z_tmp[,1]
  z <- z_tmp[,2]
  
  # dimensions
  m <- length(s.id)
  n <- length(id)
  
  z_s <- rep(NA, m)
  
  for (g in id)
    z_s[s.id == g] <- z[id == g]
  
  if (!is.null(x)) {
    
    x <- as.matrix(x)
    p <- ncol(design)
    
    mu_z <- sum(wts*z)/m
    mu_x <- colMeans(design)
    nu <- m - sum(wts^2)/m
    tau2 <- sum((s - z_s)^2)/(m - n)
    sigma2 <- (sum(wts*(z - mu_z)^2) - (n - 1)*tau2)/nu
    psi <- t(wts*(z - mu_z))%*%(design - matrix(rep(mu_x, n), nrow = n, byrow = TRUE))/nu
    Omega <- cov(design)
    
    phi <- c(sigma2, psi)
    Sigma <- matrix(NA, p + 1, p + 1)
    Sigma[2:(p+1),] <- Sigma[,2:(p+1)] <- psi 
    Sigma[2:(p+1),2:(p+1)] <- Omega
    
    a <- sapply(1:n, function(i, ...) {
      
      Sigma[1,1] <- sigma2 + tau2/wts[i]
      star <- c(z[i] - mu_z, design[i,] - mu_x)
      out <- c(mu_z + t(phi)%*%solve(Sigma)%*%star)
      
      return(out)
      
    })
    
  } else {
    
    mu_z <- sum(wts*z)/m
    nu <- m - sum(wts^2)/m
    tau2 <- sum((s - z_s)^2)/(m - n)
    sigma2 <- (sum(wts*(z - mu_z)^2) - (n - 1)*tau2)/nu
    
    a <- sapply(1:n, function(i, ...) {
      
      out <- c(mu_z + sigma2/(sigma2 + tau2/wts[i])*c(z[i] - mu_z))
      
      return(out)
      
    })
    
  }
    
  return(a)
  
}

agg <- function(x, id) {
  
  n <- nrow(x)
  p <- ncol(x)
  
  # initialize exposures
  z_tmp <- aggregate(x, by = list(id), mean)
  idx <- z_tmp[,1]
  z <- as.matrix(z_tmp[,2:(ncol(x)+1)])
  z_y <- matrix(NA, nrow = n, ncol = p)
  m <- nrow(z)
  wts <- c(unname(table(id)))
  
  for (g in idx)
    z_y[id == g,] <- z[idx == g,]
  
  mu_z <- colSums(wts*z)/n
  nu <- n - sum(wts^2)/n
  
  Omega.tmp <- matrix(0, nrow = p, ncol = p)
  Sigma.tmp <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:n)
    Omega.tmp <- Omega.tmp + tcrossprod(x[i,] - z_y[i,])
  
  for (j in 1:m)
    Sigma.tmp <- Sigma.tmp + wts[j]*tcrossprod(z[j,] - mu_z)
  
  Omega <- Omega.tmp/(n - m) 
  Sigma <- (Sigma.tmp - (m - 1)*Omega)/nu
  
  a <- t(sapply(1:m, function(j, ...) {
    
    out <- c(mu_z + Sigma %*% solve(Sigma + Omega/wts[j])%*%c(z[j,] - mu_z))
    return(out)
    
  }))
  
  colnames(a) <- colnames(x)
  agg <- data.frame(id = idx, a = a)
  
  return(agg)
  
}

pred <- function(s, star, w, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "SL.gam")){
  
  # set up evaluation points & matrices for predictions
  ws <- data.frame(w, star)
  ws.tmp <- data.frame(ws[!is.na(s),])
  s.tmp <- s[!is.na(s)]
  colnames(ws.tmp) <- colnames(ws) <- c(colnames(w), "expos")
  
  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = s.tmp, X = ws.tmp, SL.library = sl.lib)
  stilde <- c(predict(mumod, newdata = ws)$pred)
  
  return(stilde)
  
}
