simex_dr <- function(z, y, x, id, y.id, sigma, 
                     n.boot = 100, degree = 2, lambda = seq(0.1, 2.1, by = 0.25),
                     a.vals = seq(min(a), max(a), length.out = 20), 
                     span = NULL, span.seq = seq(0.5, 1, by = 0.05), k = 10,
                     sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "sL.gam")){
  
  if (any(duplicated(id)))
    stop("duplicate id detected")
  
  size <- unname(table(y.id))
  z <- z[order(id)]
  id <- id[order(id)]
  
  eps <- replicate(n.boot, rnorm(length(id), 0, sig_epe))
  
  l.vals <- mclapply.hack(lambda, function(lam, ...){
    
    z.mat <- z + sqrt(lam)*eps
    a.mat <- matrix(NA, nrow = length(y.id), ncol = n.boot)
    
    for(g in id)
      a.mat[y.id == g,] <- z.mat[id == g,]
    
    vals <- apply(a.mat, 2, hct_dr, y = y, x = x, y.id = y.id, a.vals = a.vals, 
                  span = span, span.seq = span.seq, k = k, sl.lib = sl.lib)
    
    mu.vals <- matrix(unlist(lapply(vals, function(x) x[1])), ncol = length(vals))
    sig.vals <- matrix(unlist(lapply(vals, function(x) x[2])), ncol = length(vals))
    
    m.vals <- matrix(rep(rowMeans(mu.vals), n.boot), nrow = length(a.vals), ncol = n.boot)
    s.hat <- rowSums((mu.vals - m.vals)^2)/(n.boot - 1)
    tau.hat <- rowMeans(sig.vals)
    
    return(list(estimate = rowMeans(mu.vals), variance = tau.hat - s.hat))
    
  })
  
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
  
  return(list(estimate = estimate, variance = variance, Psi = Psi, Phi = Phi, lambda = lambda))
  
}