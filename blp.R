
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
    p <- ncol(x)
    
    mu_z <- sum(wts*z)/m
    mu_x <- colMeans(x)
    nu <- m - sum(wts^2)/m
    tau2 <- sum((s - z_s)^2)/(m - n)
    sigma2 <- (sum(wts*(z - mu_z)^2) - (n - 1)*tau2)/nu
    psi <- t(wts*(z - mu_z))%*%(x - matrix(rep(mu_x, n), nrow = n, byrow = TRUE))/nu
    Omega <- cov(x)
    
    phi <- c(sigma2, psi)
    Sigma <- matrix(NA, p + 1, p + 1)
    Sigma[2:(p+1),] <- Sigma[,2:(p+1)] <- psi 
    Sigma[2:(p+1),2:(p+1)] <- Omega
    
    a <- sapply(1:n, function(i, ...) {
      
      Sigma[1,1] <- sigma2 + tau2/wts[i]
      star <- c(z[i] - mu_z, x[i,] - mu_x)
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

pred <- function(s, star, w, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "SL.gam")){
  
  # set up evaluation points & matrices for predictions
  ws <- data.frame(w, star)
  ws.tmp <- data.frame(ws[!is.na(s),])
  s.tmp <- s[!is.na(s)]
  colnames(ws.tmp) <- colnames(ws) <- c(colnames(w), "expos")
  
  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = s.tmp, X = ws.tmp, SL.library = sl.lib)
  stilde <- c(predict(mumod, newdata = ws)$pred)
  stilde[!is.na(s)] <- s[!is.na(s)]
  
  return(stilde)
  
}
