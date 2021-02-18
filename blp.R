
blp <- function(s, x, s.id, y.id, fmla) {
  
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
  size <- table(s.id)
  
  d <- aggregate(model.matrix(as.formula(fmla), data = data.frame(x)), by = list(y.id), mean)
  design <- as.matrix(d[,2:ncol(d)])[,-1]
  
  # dimensions
  l <- length(s.id)
  m <- length(id)
  n <- length(y.id)
  p <- ncol(design)
  
  # initialize exposures
  z <- aggregate(s, by = list(s.id), mean)[,2]
  z_s <- rep(NA, length(s.id))
  
  for (g in id)
    z_s[s.id == g] <- z[id == g]
  
  mu_z <- sum(size*z)/l
  mu_x <- colMeans(design)
  nu <- l - sum(size^2)/l
  tau2 <- sum((s - z_s)^2)/(l - m)
  sigma2 <- (sum(size*(z - mu_z)^2) - (m - 1)*tau2)/nu
  psi <- t(size*(z - mu_z))%*%(design - matrix(rep(mu_x, m), nrow = m, byrow = TRUE))/nu
  Omega <- cov(design)
  
  phi <- c(sigma2, psi)
  Sigma <- matrix(NA, p + 1, p + 1)
  Sigma[2:(p+1),] <- Sigma[,2:(p+1)] <- psi 
  Sigma[2:(p+1),2:(p+1)] <- Omega
  
  a <- sapply(1:m, function(i, ...) {
    
    Sigma[1,1] <- sigma2 + tau2/size[i]
    star <- c(z[i] - mu_z, design[i,] - mu_x)
    c(mu_z + t(phi)%*%solve(Sigma)%*%star)
    
    
  })
  
  return(a)
  
}
