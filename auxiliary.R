
# regression calibration function
pred <- function(s, s.tilde, w, sl.lib = c("SL.mean", "SL.glm", "SL.glmnet", "SL.earth", "SL.ranger")){
  
  # set up evaluation points & matrices for predictions
  ws <- data.frame(w, s.tilde)
  ws.obs <- data.frame(ws[!is.na(s),])
  s.obs <- s[!is.na(s)]
  colnames(ws.obs) <- colnames(ws) <- c(colnames(w), "expos")
  
  # estimate nuisance outcome model with SuperLearner
  mumod <- SuperLearner(Y = s.obs, X = ws.obs, SL.library = sl.lib)
  s.hat<- c(predict(mumod, newdata = ws)$pred)
  s.hat[!is.na(s)] <- s[!is.na(s)]
  
  return(s.hat)
  
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
