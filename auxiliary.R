
# validation regression calibration function
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

# Check the adjacency matrix for spatial gps
quants <- function(V) {
  
  #### Check V is a matrix of the correct dimension
  if(!is.matrix(V)) stop("V is not a matrix.", call.=FALSE)
  n <- nrow(V)
  if(ncol(V)!= n) stop("V is not a square matrix.", call.=FALSE)    
  
  #### Check validity of inputed V matrix
  if(sum(is.na(V))>0) stop("V has missing 'NA' values.", call.=FALSE)
  if(!is.numeric(V)) stop("V has non-numeric values.", call.=FALSE)
  if(min(V)<0) stop("V has negative elements.", call.=FALSE)
  if(sum(V!=t(V))>0) stop("V is not symmetric.", call.=FALSE)
  if(min(apply(V, 1, sum))==0) stop("V has some areas with no neighbours (one of the row sums equals zero).", call.=FALSE)    
  
  #### Create the triplet form
  ids <- which(V > 0, arr.ind = T)
  V.triplet <- cbind(ids, V[ids])
  V.triplet <- V.triplet[ ,c(2,1,3)]
  
  n.triplet <- nrow(V.triplet) 
  V.triplet.sum <- tapply(V.triplet[ ,3], V.triplet[ ,1], sum)
  n.neighbours <- tapply(V.triplet[ ,3], V.triplet[ ,1], length)
  
  #### Create the start and finish points for V updating
  V.begfin <- cbind(c(1, cumsum(n.neighbours[-n])+1), cumsum(n.neighbours))
  
  #### Return the critical quantities
  results <- list(V = V, V.triplet = V.triplet, n.triplet = n.triplet, 
                  V.triplet.sum = V.triplet.sum, n.neighbours = n.neighbours, 
                  V.begfin = V.begfin, n = n)
  return(results)
  
}
