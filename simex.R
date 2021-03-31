
# wrapper for hct-dr.R to allow for SIMEX correction
simex <- function(s, y, x, id, s.id, family = gaussian(), offset = rep(0, length(y)), # measurement error values
                     n.boot = 100, degree = 3, lambda = seq(0.1, 2.1, by = 0.25), # simex parameters
                     a.vals = seq(min(a), max(a), length.out = 20), mc.cores = 3, # erc() parameters
                     span = NULL, span.seq = seq(0.15, 1, by = 0.05), k = 5){ # cross validation shenanigans

  
  if (any(duplicated(id)))
    stop("duplicate id detected")
  
  l.vals <- mclapply.hack(lambda, function(lam, s, y, x.mat, id, s.id, ...){
    
    z.mat <- replicate(n.boot, sapply(id, function(i, ...)
      mean(s[s.id == i]) + sqrt(lam/sum(s.id == i))*contr(s[s.id == i])))
    
    vals <- apply(z.mat, 2, erc, y = y, x = x.mat, offset = offset, family = family,
                  a.vals = a.vals, sl.lib = sl.lib, span = span, span.seq = span.seq, k = k)
    
    mu.vals <- matrix(unlist(lapply(vals, function(r) r[1])), ncol = length(vals))
    sig.vals <- matrix(unlist(lapply(vals, function(r) r[2])), ncol = length(vals))
    
    m.vals <- matrix(rep(rowMeans(mu.vals), n.boot), nrow = length(a.vals), ncol = n.boot)
    s.hat <- rowSums((mu.vals - m.vals)^2)/(n.boot - 1)
    tau.hat <- rowMeans(sig.vals)
    
    return(list(estimate = rowMeans(mu.vals), variance = tau.hat - s.hat))
    
  }, s = s, y = y, x.mat = x, id = id, s.id = s.id, mc.cores = mc.cores)
  
  if (any(lambda == 0)){
    
    Psi <- t(matrix(unlist(lapply(l.vals, function(x) x$estimate)), ncol = length(l.vals)))[,-which(lambda == 0)]
    Phi <- t(matrix(unlist(lapply(l.vals, function(x) x$variance)), ncol = length(l.vals)))[,-which(lambda == 0)]
    
    
  } else {
    
    Psi <- t(matrix(unlist(lapply(l.vals, function(x) x$estimate)), ncol = length(l.vals)))
    Phi <- t(matrix(unlist(lapply(l.vals, function(x) x$variance)), ncol = length(l.vals)))
    
  }
  
  L <- cbind(1, poly(lambda, degree = degree, raw = TRUE))
  chi <- c(1, poly(-1, degree = degree, raw = TRUE))
  estimate <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Psi)
  variance <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% Phi)
  
  out <- list(estimate = estimate, variance = variance, Psi = Psi, Phi = Phi, lambda = lambda)
  
  return(out)
  
}

contr <- function(si) {
  
  xi <- rnorm(length(si), 0, 1)
  
  if (length(xi) == 1){
    
    return(0)
    
  }
  
  ci <- (xi - mean(xi))/sqrt(sum((xi - mean(xi))^2))
  return(sum(ci*si))
    
}

