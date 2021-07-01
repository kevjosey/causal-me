
blp <- function(s, s.id, x = NULL, id = NULL) {
  
  wts <- c(unname(table(s.id)))
  
  if (!is.null(id)) {
    
    if(sum(duplicated(id)) != 0)
      stop("duplicate id detected")
    
    x <- x[order(id),]
    id.check <- id[order(id)]
    
  } else 
    id.check <- NULL
  
  # initialize exposures
  z_tmp <- aggregate(s, by = list(s.id), mean)
  id <- z_tmp[,1]
  z <- z_tmp[,2]
  
  if (!is.null(id.check)) {
    
    if (!identical(id, id.check))
      stop("id provided does not match unique values of s.id")
    
    if (nrow(x) != length(id))
      stop("x provided does not match length of unique s.id")
    
  }
  
  # dimensions
  m <- length(s.id)
  n <- length(id)
  
  ord <- order(s.id)
  z_s <- rep(z, wts)[order(ord)]
  
  if (!is.null(x)) {
    
    x <- as.matrix(x)
    p <- ncol(x)
    
    mu_z <- sum(wts*z)/m
    mu_x <- colMeans(x)
    muMat_x <- matrix(rep(mu_x, n), nrow = n, byrow = TRUE)
    nu <- m - sum(wts^2)/m
    
    tau2 <- sum((s - z_s)^2)/(m - n)
    sigma2 <- c(sum(wts*(z - mu_z)^2) - (n - 1)*tau2)/nu
    psi <- c(crossprod(wts*(z - mu_z), as.matrix(x - muMat_x)))/nu
    Omega <- cov(x)
    
    phi <- c(sigma2, psi)
    V <- matrix(NA, p + 1, p + 1)
    V[2:(p+1),] <- V[,2:(p+1)] <- psi 
    V[2:(p+1),2:(p+1)] <- Omega
    
    a <- sapply(1:n, function(i, ...) {
      
      V[1,1] <- sigma2 + tau2/wts[i]
      star <- c(z[i] - mu_z, t(x[i,]) - mu_x)
      out <- c(mu_z + t(phi)%*%solve(V)%*%star)
      
      return(out)
      
    })
    
    cvar <- sapply(1:n, function(i, ...) {
      
      V[1,1] <- sigma2 + tau2/wts[i]
      out <- c(sigma2 - t(phi)%*%solve(V)%*%phi)
      
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
    
    cvar <- sapply(1:n, function(i, ...) {
      
      V[1,1] <- sigma2 + tau2/wts[i]
      out <- c(sigma2 - (sigma2)^2/(sigma2 + tau2/wts[i]))
      
      return(out)
      
    })
    
  }
  
  return(list(id = id, a = a, cvar = cvar))
  
}

multi_blp <- function(s, s.id, x = NULL, id = NULL) {
  
  wts <- c(unname(table(s.id)))
  
  if (!is.null(id)) {
    
    if(sum(duplicated(id)) != 0)
      stop("duplicate id detected")
    
    x <- x[order(id),]
    id.check <- id[order(id)]
    
  } else 
    id.check <- NULL
  
  # initialize exposures
  z_tmp <- aggregate(s, by = list(s.id), mean)
  id <- z_tmp[,1]
  z <- as.matrix(z_tmp[,2:(ncol(s)+1)])
  
  if (!is.null(id.check)) {
    
    if (!identical(id, id.check))
      stop("id provided does not match unique values of s.id")
    
    if (nrow(x) != length(id))
      stop("x provided does not match length of unique s.id")
    
  }
  
  # dimensions
  m <- length(s.id)
  n <- length(id)
  
  ord <- order(s.id)
  z_tmp <- z[rep(1:nrow(z), wts),]
  z_s <- z_tmp[order(ord),]
  
  p <- ncol(s)
  q <- ncol(x)
  
  if (!is.null(x)) {
    
    mu_z <- colSums(wts*z)/m
    mu_x <- colMeans(x)
    muMat_x <- matrix(rep(mu_x, n), nrow = n, byrow = TRUE)
    muMat_z <- matrix(rep(mu_z, n), nrow = n, byrow = TRUE)
    nu <- m - sum(wts^2)/m
    
    Omega <- crossprod(as.matrix(s - z_s))/(m - n)
    Sigma <- as.matrix(crossprod(wts*(z - muMat_z), (z - muMat_z)) - (n - 1)*Omega)/nu
    Psi <- crossprod(wts*(z - muMat_z), as.matrix(x - muMat_x))/nu
    Tau <- cov(x)
    
    Phi <- cbind(Sigma, Psi)
    V <- matrix(NA, p + q, p + q)
    V[(p+1):(p+q),1:p] <- V[1:p,(p+1):(p+q)] <- Psi
    V[(p+1):(p+q),(p+1):(p+q)] <- Tau
    
    a <- t(sapply(1:n, function(i, ...) {
      
      V[1:p,1:p] <- Sigma + Omega/wts[i]
      star <- c(t(z[i,]) - mu_z, t(x[i,]) - mu_x)
      out <- c(mu_z + Phi %*% solve(V) %*% star)
      
      return(out)
      
    }))
    
  } else {
    
    mu_z <- colSums(wts*z)/m
    muMat_z <- matrix(rep(mu_z, n), nrow = n, byrow = TRUE)
    nu <- m - sum(wts^2)/m
    
    Omega <- crossprod(as.matrix(s - z_s))/(m - n)
    Sigma <- as.matrix(crossprod(wts*(z - muMat_z), (z - muMat_z)) - (n - 1)*Omega)/nu
    
    a <- t(sapply(1:n, function(i, ...) {
      
      V <- Sigma + Omega/wts[i]
      out <- c(mu_z + Sigma %*% solve(V) %*% c(t(z[i,]) - mu_z))
      
      return(out)
      
    }))
    
  }
  
  colnames(a) <- colnames(x)
  return(data.frame(id = id, a))
  
}

pred <- function(s, star, w, sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.earth")){
  
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