## Inputs

# n = total sample size
# d = observed true sample size d < n
# sig_gps = sd of the true exposure
# sig_epe = sd of the error prone exposure
# sig_out = sd of the outcome
# epe_scen = error prone exposure misspecification
# gps_scen = gpse exposure misspecification
# out_scen = outcome misspecification

## Output

# s = observed exposure, either true or error prone
# t_obs = true exposures when they are observed
# y = outcome
# group = cluster assignment
# X = covariate matrix
# gps = generalized propensity score of the true exposure

gen_agg_data <- function(m, n, prob = 0.1) {
  
  if (m > n)
    stop("you stop that, you!")
  
  # covariates
  u1 <- stats::rnorm(m, 0, 1)
  u2 <- stats::rnorm(m, 0, 1)
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- x4 <- rep(NA, n)
  
  y.id <- sample(1:m, n, replace = TRUE)
  s.id <- 1:m
  
  for (g in s.id) {
    
    x3[y.id == g] <- u1[s.id == g]
    x4[y.id == g] <- u2[s.id == g]
    
  }
  
  x <- cbind(x1, x2, x3, x4)
  v <- aggregate(x, by = list(y.id), mean)[,2:(ncol(x) + 1)]
  
  mu_gps <- 1 + 0.5*v[,1] - 0.5*v[,2] - 0.5*v[,3] + 0.5*v[,4]
  
  a <- rnorm(m, mu_gps, sig_gps)
  z <- rnorm(m, a, sig_epe)
  a_y <- rep(NA, n)
  z_y <- rep(NA, n)
  
  for (g in s.id){
    a_y[y.id == g] <- a[s.id == g]
    z_y[y.id == g] <- z[s.id == g]
  }
  
  mu_out <- -0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4 + z_y*(1 + 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  
  y <- rbinom(n, 1, plogis(mu_out))
  
  gps <- dnorm(a, mu_gps, sig_gps)
  
  # create simulation dataset
  sim <- list(a = a, z = z, y = y, x = x, s.id = s.id, y.id = y.id, gps = gps, a_y = a_y, z_y = z_y)
  
  return(sim)
  
}

predict.example <- function(a, x, id){
  
  v3 <- unique(x[,3])[order(unique(id))]
  v4 <- unique(x[,4])[order(unique(id))]
  
  # transformed predictors
  u1 <- as.numeric(scale((x[,1] + x[,2])^2))
  u2 <- as.numeric(scale(cos(2*x[,2])))
  z3 <- as.numeric(scale(sin(2*v3)))
  z4 <- as.numeric(scale(-abs(v3 + v4)))
  
  lvl <- unique(id)[order(unique(id))]
  u3 <- u4 <- rep(NA, length(id))
  
  for (g in lvl) {
    
    u3[id == g] <- z3[lvl == g]
    u4[id == g] <- z4[lvl == g]
    
  }
  
  u <- cbind(u1, u2, u3, u4)
  
  out <- rep(NA, length(a))
  gamma <- unname(table(y.id)/mean(table(y.id)))
  
  for(i in 1:length(a)) {
    
    mu_out <- aggregate(x %*% c(-0.75,-0.25,0.25,0.75) + rep(a[i],length(y.id))*(1 + x %*% c(0.5, -0.5, 0.5, -0.5)), 
                        by = list(id), mean)[,2]
    out[i] <- mean(gamma*plogis(mu_out))
    
  }
  
  return(out)
  
}
