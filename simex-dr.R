simex_dr <- function(z, y, x, s.id, y.id, sig_epe, n.boot = 100, degree = 2,
                     lambda = seq(0, 5, by = 0.25),
                     a.vals = seq(min(a), max(a), length.out = 20), 
                     span = NULL, span.seq = seq(0.5, 1, by = 0.05), k = 10,
                     sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth", "sL.gam")){
  
  # first some house-cleaning
  lvl_y <- unique(y.id)[order(unique(y.id))]
  lvl_s <- s.id[order(s.id)]
  z <- z[order(s.id)]
  
  if(length(lvl_s) != length(lvl_y))
    stop("some levels in y.id are not represented by s.id. There is no exposure data for these entries.")
  
  if(!all(lvl_s == lvl_y))
    stop("some levels in y.id are not represented by s.id. There is no exposure data for these entries.")
  
  z <- z[lvl_s %in% lvl_y]
  s.id <- lvl_s[lvl_s %in% lvl_y]
  
  eps <- replicate(n.boot, rnorm(m, 0, sig_epe))
  
  l.vals <- sapply(lambda, function(lam, ...){
    
    z.mat <- z + sqrt(lam)*eps
    a.mat <- matrix(NA, nrow = n, ncol = n.boot)
    
    for(g in lvl_y)
      a.mat[y.id == g,] <- z.mat[lvl_y == g,]
    
    vals <- apply(a.mat, 2, hct_dr, y = y, x = x, y.id = y.id, a.vals = a.vals, 
                  span = span, span.seq = span.seq, k = k, sl.lib = sl.lib)
    
    mu.vals <- matrix(unlist(lapply(vals, function(x) x[1])), ncol = length(vals))
    sig.vals <- matrix(unlist(lapply(vals, function(x) x[2])), ncol = length(vals))
    
    m.vals <- matrix(rep(rowMeans(mu.vals), n.boot), nrow = length(a.vals), ncol = n.boot)
    s.hat <- rowSums((mu.vals - m.vals)^2)/(n.boot - 1)
    tau.hat <- rowMeans(sig.vals)
    
    return(list(estimate = m.vals[,1], variance = tau.hat - s.hat))
    
  })
  
  if (any(lambda == 0)){
    
    Psi <- do.call(cbind, l.vals[1,])[,-which(lambda == 0)]
    Phi <- do.call(cbind, l.vals[2,])[,-which(lambda == 0)]
    
  } else {
    
    Psi <- do.call(cbind, l.vals[1,])
    Phi <- do.call(cbind, l.vals[2,])
    
  }
  
  L <- cbind(1, poly(lambda, degree = degree, raw = TRUE))
  chi <- c(1, poly(-1, degree = degree, raw = TRUE))
  estimate <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% t(Psi))
  variance <- c(t(chi) %*% solve(t(L) %*% L) %*% t(L) %*% t(Phi))
  
  return(list(estimate = estimate, variance = variance, Psi = Psi, Phi = Phi, lambda = lambda))
  
}