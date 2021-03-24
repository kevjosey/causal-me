## Inputs

# n = number of clusters
# m = number of grids m < n
# l = number of individuals comprising clusters
# sig_gps = sd of the true exposure
# sig_agg = sd of the error prone exposure
# gps_scen = gps exposure misspecification
# out_scen = outcome misspecification

## Output

# s = error prone exposure (on a grid)
# a = true exposure (on cluster)
# a_y = true cluster exposure expanded to individuals
# y = outcome (individual)
# y.id = cluster assignment for y
# s.id = cluster assignment for s
# id = id for a corresponding to y.id
# x = covariate matrix

gen_data <- function(m, n, sig_agg = sqrt(2), sig_gps = 1, sig_pred = sqrt(0.5),
                     gps_scen = c("a", "b"), out_scen = c("a", "b"), pred_scen = c("a", "b")) {
  
  if (n > m)
    stop("you stop that, you!")
    
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  
  w1 <- stats::rnorm(m, 1, 2)
  w2 <- rep(NA, m)
  
  id <- 1:n
  s.id <- sample(id, m, replace = TRUE)
  offset <- runif(n, 5, 50)
  
  # transformed predictors
  u1 <- as.numeric(scale((x1 + x2)^2))
  u2 <- as.numeric(scale(cos(2*x2)))
  u3 <- as.numeric(scale(sin(2*x3)))
  u4 <- as.numeric(scale(-abs(x3 + x4)))
  
  if (gps_scen == "b") {
    
    for (g in 1:n)
      w2[s.id == g] <- u4[g] + rnorm(sum(s.id == g), 0, 1)
    
    mu_gps <- 8 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4
    
    
  } else {
    
    for (g in 1:n)
      w2[s.id == g] <- x4[g] + rnorm(sum(s.id == g), 0, 1)
    
    mu_gps <- 8 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4
    
  }
  
  x <- cbind(x1, x2, x3, x4)
  u <- cbind(u1, u2, u3, u4)
  w <- cbind(w1, w2)
  
  a <- rnorm(n, mu_gps, sig_gps)
  a_s <- rep(NA, m)
  
  for (g in 1:n)
    a_s[s.id == g] <- a[g]
    
  s <- rnorm(m, a_s, sig_agg)
  
  if (pred_scen == "b"){
    star <- rnorm(m, s - 0.1 + 0.5*w1 - 0.5*w2, sig_pred)
  } else {
    star <- rnorm(m, s, sig_pred)
  }
  
  if (out_scen == "b") {
    mu_out <- -3 - 0.3*u1 - 0.1*u2 + 0.1*u3 + 0.3*u4 + 0.3*(a - 8) - 0.1*(a - 8)^2
  } else { # y_scen == "b"
    mu_out <- -3 - 0.3*x1 - 0.1*x2 + 0.1*x3 + 0.3*x4 + 0.3*(a - 8) - 0.1*(a - 8)^2
  }
  
  y <- rpois(n, exp(mu_out + log(offset)))
  
  # create simulation dataset
  sim <- list(a = a, s = s, star = star, y = y, x = x, w = w, 
              id = id, s.id = s.id, offset = offset)
  return(sim)
  
}

predict_example <- function(a.vals, x, out_scen = c("a", "b")) {
  
  # transformed predictors
  u1 <- as.numeric(scale((x[,1] + x[,1])^2))
  u2 <- as.numeric(scale(cos(2*x[,2])))
  u3 <- as.numeric(scale(sin(2*x[,3])))
  u4 <- as.numeric(scale(-abs(x[,3] + x[,4])))
  
  u <- cbind(u1, u2, u3, u4)
  out <- rep(NA, length(a.vals))
  
  for(i in 1:length(a.vals)) {
    
    a.vec <- rep(a.vals[i],nrow(x))
    
    if (out_scen == "b") {
      mu_out <- exp(-3 + u %*% c(-0.3,-0.1,0.1,0.3) + 0.3*(a.vec - 8) - 0.1*(a.vec - 8)^2)
    } else { # out_scen == "a"
      mu_out <- exp(-3 + x %*% c(-0.3,-0.1,0.1,0.3) + 0.3*(a.vec - 8) - 0.1*(a.vec - 8)^2)
    }
    
    out[i] <- mean(mu_out)
    
  }
  
  return(out)
  
}
