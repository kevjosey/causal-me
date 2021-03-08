
blp <- function(s, s.id, y.id, x = NULL) {
  
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
  wts <- table(s.id)
  
  if (!is.null(x)) {
    
    d <- aggregate(model.matrix(~ ., data = data.frame(x)), by = list(y.id), mean)
    design <- as.matrix(d[,2:ncol(d)])[,-1]
    p <- ncol(design)
    
  }
  
  # dimensions
  l <- length(s.id)
  m <- length(id)
  n <- length(y.id)

  # initialize exposures
  z <- aggregate(s, by = list(s.id), mean)[,2]
  z_s <- rep(NA, length(s.id))
  
  for (g in id)
    z_s[s.id == g] <- z[id == g]
  
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
    
  return(list(a = a, a_y = a_y, a_s = a_s))
  
}
