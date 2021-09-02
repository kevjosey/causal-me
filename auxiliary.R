
ipw <- function(a, x, beta, sigma2, a.vals) {
  
  n <- length(a)
  x.new <- x[rep(1:n, length(a.vals) + 1), ]
  a.new <- c(a, rep(a.vals, each = n))
  pimod.vals <- c(x.new %*% beta)
  pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2))
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
  phat[which(phat < 0)] <- 1e-6
  out <- phat/pihat
  return(out)
  
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

# highest posterior density
hpd <- function(x, alpha = 0.05){
  
  n <- length(x)
  m <- round(n * alpha)
  x <- sort(x)
  y <- x[(n - m + 1):n] - x[1:m]
  z <- min(y)
  k <- which(y == z)[1]
  c(x[k], x[n - m + k])
  
}
