
blp <- function(s, s.id, x = NULL) {
  
  wts <- c(unname(table(s.id)))
  
  # initialize exposures
  z_tmp <- aggregate(s, by = list(s.id), mean)
  id <- z_tmp[,1]
  z <- z_tmp[,2]
  
  # dimensions
  l <- length(s.id)
  m <- length(id)
  
  z_s <- rep(NA, length(s.id))
  
  for (g in id)
    z_s[s.id == g] <- z[id == g]
  
  if (!is.null(x)) {
    
    d <- aggregate(model.matrix(~ ., data = data.frame(x)), by = list(y.id), mean)
    design <- as.matrix(d[,2:ncol(d)])[,-1]
    p <- ncol(design)
    
  }
  
  if (!is.null(x)) {
    
    mu_z <- sum(wts*z)/l
    mu_x <- colMeans(design)
    nu <- l - sum(wts^2)/l
    tau2 <- sum((s - z_s)^2)/(l - m)
    sigma2 <- (sum(wts*(z - mu_z)^2) - (m - 1)*tau2)/nu
    psi <- t(wts*(z - mu_z))%*%(design - matrix(rep(mu_x, m), nrow = m, byrow = TRUE))/nu
    Omega <- cov(design)
    
    phi <- c(sigma2, psi)
    Sigma <- matrix(NA, p + 1, p + 1)
    Sigma[2:(p+1),] <- Sigma[,2:(p+1)] <- psi 
    Sigma[2:(p+1),2:(p+1)] <- Omega
    
    a <- sapply(1:m, function(i, ...) {
      
      Sigma[1,1] <- sigma2 + tau2/wts[i]
      star <- c(z[i] - mu_z, design[i,] - mu_x)
      out <- c(mu_z + t(phi)%*%solve(Sigma)%*%star)
      
      return(out)
      
    })
    
  } else {
    
    mu_z <- sum(wts*z)/l
    nu <- l - sum(wts^2)/l
    tau2 <- sum((s - z_s)^2)/(l - m)
    sigma2 <- (sum(wts*(z - mu_z)^2) - (m - 1)*tau2)/nu
    
    a <- sapply(1:m, function(i, ...) {
      
      out <- c(mu_z + sigma2/(sigma2 + tau2/wts[i])*c(z[i] - mu_z))
      
      return(out)
      
    })
    
  }
  
  a_s <- rep(NA, length(s.id))
  a_y <- rep(NA, length(y.id))
  
  for (g in id) {
    
    a_s[s.id == g] <- a[id == g]
    a_y[y.id == g] <- a[id == g]
    
  }
    
  return(list(a = a, sigma2 = sigma2, a_s = a_s))
  
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

pred <- function(s, shat, w, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "SL.gam")){
  
  # set up evaluation points & matrices for predictions
  ws <- data.frame(w, shat)
  ws.tmp <- data.frame(ws[!is.na(s),])
  s.tmp <- s[!is.na(s)]
  colnames(ws.tmp) <- colnames(ws) <- c(colnames(w), "expos")
  
  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = s.tmp, X = ws.tmp, SL.library = sl.lib)
  stilde <- c(predict(mumod, newdata = ws)$pred)
  
  return(stilde)
  
}
