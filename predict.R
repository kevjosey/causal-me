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

  
