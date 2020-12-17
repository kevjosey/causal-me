gen_agg_dr <- function(n, m, sig_gps, gps_scen = c("a", "b"), out_scen = c("a", "b")) {
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  v3 <- stats::rnorm(m, 0, 1)
  v4 <- stats::rnorm(m, 0, 1)
  
  x3 <- x4 <- rep(NA, n)
  
  id <- sample(1:m, n, replace = TRUE)
  lvl <- unique(id)[order(unique(id))]
  
  for (g in lvl) {
    
    x3[id == g] <- v3[lvl == g]
    x4[id == g] <- v4[lvl == g]
    
  }
  
  # transformed predictors
  u1 <- as.numeric(scale((x1 + x2)^2))
  u2 <- as.numeric(scale(cos(2*x2)))
  z3 <- as.numeric(scale(sin(2*v3)))
  z4 <- as.numeric(scale(-abs(v3 + v4)))
  
  u3 <- u4 <- rep(NA, n)
  
  for (g in lvl) {
    
    u3[id == g] <- z3[lvl == g]
    u4[id == g] <- z4[lvl == g]
    
  }
  
  # propensity score
  if (gps_scen == "b") {
    mu_gps <- 1 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4
  } else { # z_scen == "a"
    mu_gps <- 1 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4
  }
  
  tg <- aggregate(rnorm(n, mu_gps, sig_gps), by = list(id), mean)[,2]
  t <- rep(NA, n)
  
  for (g in lvl)
    t[id == g] <- tg[lvl == g]
  
  if (out_scen == "b") {
    mu_out <- -0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4 + t*(1 + 0.5*u1 - 0.5*u2 + 0.5*u3 - 0.5*u4)
  } else { # y_scen == "b"
    mu_out <- -0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4 + t*(1 + 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  }
  
  y <- rbinom(n, 1, plogis(mu_out))
  
  X <- cbind(x1, x2, x3, x4)
  gps <- aggregate(dnorm(t, mu_gps, sig_gps), by = list(id), mean)[,2]
  
  # create simulation dataset
  sim <- list(y = y, t = t, X = X, id = id, gps = gps)
  
  return(sim)
  
}

predict.example <- function(a, x, id, out_scen = c("a", "b")){
  
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
  
  for(i in 1:length(a)) {
    
    if (out_scen == "b") {
      mu_out <- u %*% c(-0.75,-0.25,0.25,0.75) + rep(a[i],n)*(1 + u %*% c(0.5, -0.5, 0.5, -0.5))
    } else { # y_scen == "b"
      mu_out <- x %*% c(-0.75,-0.25,0.25,0.75) + rep(a[i],n)*(1 + x %*% c(0.5, -0.5, 0.5, -0.5))
    }
    
    out[i] <- mean(plogis(mu_out))
    
  }
  
  return(out)
  
}
