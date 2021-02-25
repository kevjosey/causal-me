
blp <- function(s, s.id, y.id, t = NULL, x = NULL, w = NULL,  
                sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "sL.gam")) {
  
  if (!is.null(w) & !is.null(t)) {
    
    # set up evaluation points & matrices for predictions
    ws <- data.frame(w, s)
    ws.tmp <- ws[!is.na(t),]
    t.tmp <- ws[!is.na(t)]
    
    # estimate nuisance outcome model with SuperLearner
    mumod <- SuperLearner(Y = t.tmp, X = ws.tmp, SL.library = sl.lib)
    s <- predict(mumod, newdata = ws)
    
  } else if (!is.null(w) & is.null(t)) {
    
    stop("!is.null(w) & is.null(t)")
    
  } else if (is.null(w) & !is.null(t)) {
    
    s[!is.na(t)] <- t[!is.na(t)]
    warning("replaced values of s with t wherever available.")
    
  }
  
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
  offset <- table(s.id)
  
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
    
    mu_z <- sum(offset*z)/l
    mu_x <- colMeans(design)
    nu <- l - sum(offset^2)/l
    tau2 <- sum((s - z_s)^2)/(l - m)
    sigma2 <- (sum(offset*(z - mu_z)^2) - (m - 1)*tau2)/nu
    psi <- t(offset*(z - mu_z))%*%(design - matrix(rep(mu_x, m), nrow = m, byrow = TRUE))/nu
    Omega <- cov(design)
    
    phi <- c(sigma2, psi)
    Sigma <- matrix(NA, p + 1, p + 1)
    Sigma[2:(p+1),] <- Sigma[,2:(p+1)] <- psi 
    Sigma[2:(p+1),2:(p+1)] <- Omega
    
    a <- sapply(1:m, function(i, ...) {
      
      Sigma[1,1] <- sigma2 + tau2/offset[i]
      star <- c(z[i] - mu_z, design[i,] - mu_x)
      c(mu_z + t(phi)%*%solve(Sigma)%*%star)
      
      
    })
    
  } else {
    
    mu_z <- sum(offset*z)/l
    nu <- l - sum(offset^2)/l
    tau2 <- sum((s - z_s)^2)/(l - m)
    sigma2 <- (sum(offset*(z - mu_z)^2) - (m - 1)*tau2)/nu
    
    a <- sapply(1:m, function(i, ...) {
      
      c(mu_z + sigma2/(sigma2 + tau2/offset[i])*c(z[i] - mu_z))
      
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
