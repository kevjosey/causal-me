## Inputs

# n = total sample size
# m = number of clusters m < n
# l = total number of observations l < n
# sig_gps = sd of the true exposure
# sig_epe = sd of the error prone exposure
# gps_scen = gpse exposure misspecification
# out_scen = outcome misspecification

## Output

# z (or s) = observed exposure, either true or error prone
# a (or t) = true exposures when its observed
# y = outcome
# y.id = cluster assignment for y
# z.id (or s.id) = cluster assignment for z
# x = covariate matrix
# gps = generalized propensity score of the true exposure

gen_simex_data <- function(m, n, sig_epe = 2, sig_gps = 1, prob = 0.1) {
  
  if (m > n)
    stop("you stop that, you!")
  
  # covariates
  u1 <- stats::rnorm(m, 0, 1)
  u2 <- stats::rnorm(m, 0, 1)
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- x4 <- rep(NA, n)
  
  y.id <- sample(1:m, n, replace = TRUE)
  z.id <- 1:m
  gamma <- unname(table(y.id))
  
  for (g in 1:m) {
    
    x3[y.id == g] <- u1[g]
    x4[y.id == g] <- u2[g]
    
  }
  
  mu_gps <- 1 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4
  a <- aggregate(rnorm(n, mu_gps, sig_gps), by = list(y.id), mean)[,2]
  z <- rnorm(m, a, sig_epe/sqrt(unname(table(y.id))))
  a_y <- rep(NA, n)
  z_y <- rep(NA, n)
  
  for (g in 1:m){
    a_y[y.id == g] <- a[g]
    z_y[y.id == g] <- z[g]
  }
  
  mu_out <- -0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4 + a_y*(1 + 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  y <- rbinom(n, 1, plogis(mu_out))
  x <- cbind(x1, x2, x3, x4)
  
  # create simulation dataset
  sim <- list(a = a, z = z, y = y, x = x, z.id = z.id, y.id = y.id, a_y = a_y, z_y = z_y)
  
  return(sim)
  
}

gen_bayes_data <- function(l, m, n, sig_epe = 2, sig_gps = 1, prob = 0.1) {
  
  if (m > n)
    stop("you stop that, you!")
  
  # covariates
  u1 <- stats::rnorm(m, 0, 1)
  u2 <- stats::rnorm(m, 0, 1)
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  
  y.id <- sample(1:m, n, replace = TRUE)
  s.id <- sample(1:m, l, replace = TRUE)
  
  x3 <- x4 <- rep(NA, n)
  
  for (g in 1:m) {
    
    x3[y.id == g] <- u1[g]
    x4[y.id == g] <- u2[g]
    
  }
  
  x <- cbind(x1, x2, x3, x4)
  v <- aggregate(x, by = list(y.id), mean)[,2:(ncol(x) + 1)]
  mu_gps <- 1 + 0.5*v[,1] - 0.5*v[,2] - 0.5*v[,3] + 0.5*v[,4]
  a <- rnorm(m, mu_gps, sig_gps)
  a_y <- rep(NA, n)
  a_s <- rep(NA, l)
  
  for (g in 1:m) {
    
    a_y[y.id == g] <- a[g]
    a_s[s.id == g] <- a[g]
    
  }
  
  s <- rnorm(l, a_s, sig_epe)

  mu_out <- -0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4 + a_y*(1 + 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  y <- rbinom(n, 1, plogis(mu_out))
  x <- cbind(x1, x2, x3, x4)
  
  # create simulation dataset
  sim <- list(s = s, t = t, y = y, x = x, s.id = s.id, y.id = y.id, a = a, a_y = a_y)
  
  return(sim)
  
}

gen_dr_data <- function(n, m, sig_gps, gps_scen = c("a", "b"), out_scen = c("a", "b")) {
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  v3 <- stats::rnorm(m, 0, 1)
  v4 <- stats::rnorm(m, 0, 1)
  
  x3 <- x4 <- rep(NA, n)
  
  id <- sample(1:m, n, replace = TRUE)
  
  for (g in 1:m) {
    
    x3[id == g] <- v3[g]
    x4[id == g] <- v4[g]
    
  }
  
  # transformed predictors
  u1 <- as.numeric(scale((x1 + x2)^2))
  u2 <- as.numeric(scale(cos(2*x2)))
  z3 <- as.numeric(scale(sin(2*v3)))
  z4 <- as.numeric(scale(-abs(v3 + v4)))
  
  u3 <- u4 <- rep(NA, n)
  
  for (g in 1:m) {
    
    u3[id == g] <- z3[g]
    u4[id == g] <- z4[g]
    
  }
  
  # propensity score
  if (gps_scen == "b") {
    mu_gps <- 1 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4
  } else { # z_scen == "a"
    mu_gps <- 1 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4
  }
  
  tg <- aggregate(rnorm(n, mu_gps, sig_gps), by = list(id), mean)[,2]
  t <- rep(NA, n)
  
  for (g in 1:m)
    t[id == g] <- tg[g]
  
  if (out_scen == "b") {
    mu_out <- -0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4 + t*(1 + 0.5*u1 - 0.5*u2 + 0.5*u3 - 0.5*u4)
  } else { # y_scen == "b"
    mu_out <- -0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4 + t*(1 + 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  }
  
  y <- rbinom(n, 1, plogis(mu_out))
  x <- cbind(x1, x2, x3, x4)
  
  # create simulation dataset
  sim <- list(y = y, t = t, x = x, id = id)
  
  return(sim)
  
}

predict_example <- function(a, x, id, out_scen = c("a", "b")){
  
  v3 <- unique(x[,3])[order(unique(id))]
  v4 <- unique(x[,4])[order(unique(id))]
  
  # transformed predictors
  u1 <- as.numeric(scale((x[,1] + x[,2])^2))
  u2 <- as.numeric(scale(cos(2*x[,2])))
  z3 <- as.numeric(scale(sin(2*v3)))
  z4 <- as.numeric(scale(-abs(v3 + v4)))
  
  lvl <- unique(id)[order(unique(id))]
  u3 <- u4 <- rep(NA, n)
  
  for (g in lvl) {
    
    u3[id == g] <- z3[lvl == g]
    u4[id == g] <- z4[lvl == g]
    
  }
  
  u <- cbind(u1, u2, u3, u4)
  
  out <- rep(NA, length(a))
  gamma <- unname(table(y.id))
  
  for(i in 1:length(a)) {
    
    if (out_scen == "b") {
      mu_out <- aggregate(plogis(u %*% c(-0.75,-0.25,0.25,0.75) + 
                                   rep(a[i],nrow(x))*(1 + u %*% c(0.5, -0.5, 0.5, -0.5))),
                          by = list(y.id), mean)[,2]
    } else { # y_scen == "b"
      mu_out <- aggregate(plogis(x %*% c(-0.75,-0.25,0.25,0.75) + 
                                   rep(a[i],nrow(x))*(1 + x %*% c(0.5, -0.5, 0.5, -0.5))),
                          by = list(y.id), mean)[,2]
    }
    
    out[i] <- sum(gamma*mu_out)/sum(gamma)
    
  }
  
  return(out)
  
}