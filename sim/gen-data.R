## Inputs

# n = number of clusters
# mult = number of grids equals n*mult
# sig_gps = sd of the true exposure
# sig_agg = sd of the aggregation error
# sig_pred = sd of the prediction error
# gps_scen = gps exposure misspecification
# out_scen = outcome misspecification
# out_scen = exposure model misspecification

## Output

# s = true grid exposure
# star = error prone exposure
# a = true cluster exposure
# y = cluster outcomes
# id = cluster assignment for y
# s.id = cluster assignment for s
# x = confounders
# w = grid-level predictors
# offset = offset

gen_data <- function(n = c(400, 800), mult = c(5, 10), sig_agg = sqrt(2), sig_gps = 1, sig_pred = sqrt(0.5),
                     gps_scen = c("a", "b"), out_scen = c("a", "b"), pred_scen = c("a", "b")) {
    
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  x <- cbind(x1, x2, x3, x4)
  # x <- x%*%solve(chol(cov(x)))
  
  id <- 1:n
  offset <- floor(runif(n, 10, 1000))
  
  if (mult == 10) {
    s.id <- rep(id, rep(c(2,4,6,8,12,14,16,18), each = n/8))
  } else if (mult == 5){
    s.id <- rep(id, rep(c(1,2,3,4,6,7,8,9), each = n/8))
  }
  
  # transformed predictors
  u1 <- as.numeric(scale(exp(x[,1]/2)))
  u2 <- as.numeric(scale(x[,2]/(1 + exp(x[,1])) + 10))
  u3 <- as.numeric(scale((x[,1]*x[,3]/25 + 0.6)^3))
  u4 <- as.numeric(scale((x[,2] + x[,4] + 20)^2))
  
  u <- cbind(u1, u2, u3, u4)
  w2 <- rep(NA, mult*n)
  
  if (gps_scen == "b") {
    
    mu_gps <- 10 + 0.5*u[,1] - 0.5*u[,2] - 0.5*u[,3] + 0.5*u[,4]
    for (g in 1:n)
      w2[s.id == g] <- rnorm(sum(s.id == g), u[id == g,2], 1)
    
  } else {
    
    mu_gps <- 10 + 0.5*x[,1] - 0.5*x[,2] - 0.5*x[,3] + 0.5*x[,4]
    for (g in 1:n)
      w2[s.id == g] <- rnorm(sum(s.id == g), x[id == g,2], 1)
    
  }

  w1 <- stats::rnorm(mult*n, 0, 1)
  
  w <- cbind(w1, w2)
  a <- rnorm(n, mu_gps, sig_gps)
  
  stab <- table(s.id)
  a_s <- rep(a, stab)  
  
  s <- rnorm(mult*n, a_s, sig_agg)
  
  if (pred_scen == "b"){
    star <- rnorm(mult*n, s - 1 + 0.5*w1 + 0.5*w2, sig_pred)
  } else {
    star <- rnorm(mult*n, s, sig_pred)
  }
  
  if (out_scen == "b") {
    mu_out <- -4 - 0.5*u[,1] - 0.25*u[,2] + 0.25*u[,3] + 0.5*u[,4] +
      0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 0.25*(a - 10)*u[,1]
  } else { # y_scen == "b"
    mu_out <- -4 - 0.5*x[,1] - 0.25*x[,2] + 0.25*x[,3] + 0.5*x[,4] +
      0.25*(a - 10) - 0.75*cos(pi*(a - 6)/4) - 0.25*(a - 10)*x[,1]
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
      mu_out <- exp(-4 + u %*% c(-0.5,-0.25,0.25,0.5) + 0.25*(a.vec - 10) - 0.75*cos(pi*(a.vec - 6)/4) - 0.25*(a.vec - 10)*u[,1])
    } else { # out_scen == "a"
      mu_out <- exp(-4 + x %*% c(-0.5,-0.25,0.25,0.5) + 0.25*(a.vec - 10) - 0.75*cos(pi*(a.vec - 6)/4) - 0.25*(a.vec - 10)*x[,1])
    }

    out[i] <- mean(mu_out)

  }
  
  # out <- exp(-3 + 0.625/2 + 0.5*(a.vals - 8) - 0.25*(a.vals - 8)^2 )
  
  return(out)
  
}

