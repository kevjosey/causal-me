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

gen_data <- function(n = c(400, 800), mult = c(5, 10), sig_agg = sqrt(2), sig_gps = 1, sig_pred = sqrt(0.5),
                     gps_scen = c("a", "b"), out_scen = c("a", "b"), pred_scen = c("a", "b")) {
    
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  x <- cbind(x1, x2, x3, x4)
  x <- x%*%solve(chol(cov(x)))
  
  id <- 1:n
  offset <- runif(n, 10, 1000)
  
  if (mult == 10) {
    s.id <- rep(id, rep(c(1,3,8,13,25), each = n/5))
  } else if (mult == 5){
    s.id <- rep(id, rep(c(1,2,3,7,12), each = n/5))
  }

  # transformed predictors
  u1 <- as.numeric(scale(exp(x[,1]/2)))
  u2 <- as.numeric(scale(x[,2]/(1 + exp(x[,1])) + 10))
  u3 <- as.numeric(scale((x[,1]*x[,3]/25 + 0.6)^3))
  u4 <- as.numeric(scale((x[,2] + x[,4] + 20)^2))
  
  u <- cbind(u1, u2, u3, u4)
  u <- u%*%solve(chol(cov(u)))
  w2 <- rep(NA, mult*n)
  
  if (gps_scen == "b") {
    
    mu_gps <- 8 + 0.5*u[,1] - 0.5*u[,2] - 0.5*u[,3] + 0.5*u[,4]
    for (g in 1:n)
      w2[s.id == g] <- rnorm(sum(s.id == g), u[id == g,2], 1)
    
  } else {
    
    mu_gps <- 8 + 0.5*x[,1] - 0.5*x[,2] - 0.5*x[,3] + 0.5*x[,4]
    for (g in 1:n)
      w2[s.id == g] <- rnorm(sum(s.id == g), x2[id == g], 1)
    
  }

  w1 <- stats::rnorm(mult*n, 1, 1)
  w <- cbind(w1, w2)
  a <- rnorm(n, mu_gps, sig_gps)
  
  stab <- table(s.id)
  a_s <- rep(a, stab)  
  
  s <- rnorm(mult*n, a_s, sig_agg)
  
  if (pred_scen == "b"){
    star <- rnorm(mult*n, s - 1 + 0.5*w1 - 0.5*w2, sig_pred)
  } else {
    star <- rnorm(mult*n, s, sig_pred)
  }
  
  if (out_scen == "b") {
    mu_out <- -3 - 0.3*u[,1] - 0.1*u[,2] + 0.1*u[,3] + 0.3*u[,4] +
      0.3*(a - 8) - 0.1*(a - 8)^2
  } else { # y_scen == "b"
    mu_out <- -3 - 0.3*x[,1] - 0.1*x[,2] + 0.1*x[,3] + 0.3*x[,4] +
      0.3*(a - 8) - 0.1*(a - 8)^2
  }
  
  y <- rpois(n, exp(mu_out + log(offset)))
  
  # create simulation dataset
  sim <- list(a = a, s = s, star = star, y = y, x = x, w = w, 
              id = id, s.id = s.id, offset = offset)
  return(sim)
  
}

predict_example <- function(a.vals, x, out_scen = c("a", "b")) {
  
  # transformed predictors
  u1 <- as.numeric(scale(exp(x[,1]/2)))
  u2 <- as.numeric(scale(x[,2]/(1 + exp(x[,1])) + 10))
  u3 <- as.numeric(scale((x[,1]*x[,3]/25 + 0.6)^3))
  u4 <- as.numeric(scale((x[,2] + x[,4] + 20)^2))

  u <- cbind(u1, u2, u3, u4)
  u <- u%*%solve(chol(cov(u)))
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
  
  # out <- exp(-3 + 0.2/2 + 0.3*(a.vals - 8) - 0.1*(a.vals - 8)^2 )
  
  return(out)
  
}
