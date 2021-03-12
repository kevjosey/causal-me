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

gen_data <- function(l, m, n, sig_agg = sqrt(2), sig_gps = 1, sig_pred = sqrt(0.5),
                     gps_scen = c("a", "b"), out_scen = c("a", "b")) {
  
  if (m > n | l > n | m > l)
    stop("you stop that, you!")
    
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  w1 <- stats::rnorm(m, 0, 1)
  
  s.id <- sample(1:n, m, replace = TRUE)
  y.id <- sample(1:n, l, replace = TRUE)
  y.id <- rep(NA, n)
  
  for (g in 1:n)
    w2[s.id == g] <- x4[g] + rnorm(sum(s.id == g), 0, 1)
  
  # transformed predictors
  u1 <- as.numeric(scale((x1 + x4)^2))
  u2 <- as.numeric(scale(cos(2*x2)))
  u3 <- as.numeric(scale(sin(2*x4)))
  u4 <- as.numeric(scale(-abs(x2 + x3)))
  
  x <- cbind(x1, x2, x3, x4)
  u <- cbind(u1, u2, u3, u4)
  w <- cbind(w1, w2)
  
  mu_gps <- 1 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4
  
  a <- rnorm(n, mu_gps, sig_gps)
  a_s <- rep(NA, m)
  
  for (g in 1:n)
    a_s[s.id == g] <- a[g]
    
  s <- rnorm(l, a_s, sig_agg)
  t <- rnorm(n, s_y, sig_berk)
  
  mu_pred <- s + 0.25*s*w1 - 0.25*s*w2
  star <- rnorm(l, mu_pred, sig_pred)
  
  if (out_scen == "b") {
    mu_out <- -5 -0.75*u1 - 0.25*u2 - 0.25*u3 + 0.75*u4 + a*(1 - 0.5*u1 - 0.5*u2 + 0.5*u3 - 0.5*u4)
  } else { # y_scen == "b"
    mu_out <- -5 -0.75*x1 - 0.25*x2 - 0.25*x3 + 0.75*x4 + a*(1 - 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  }
  
  y <- rpois(n, exp(mu_out))
  
  # create simulation dataset
  sim <- list(a = a, s = s, star = star, y = y, x = x, w = w, s.id = s.id, )
  
  return(sim)
  
}

gen_dr_data <- function(n, m, sig_gps = 1, gps_scen = c("a", "b"), out_scen = c("a", "b")) {
  
  if (n > m)
    stop("you stop that, you!")
  
  # covariates
  x1 <- stats::rnorm(m, 0, 1)
  x2 <- stats::rnorm(m, 0, 1)
  v3 <- stats::rnorm(n, 0, 1)
  v4 <- stats::rnorm(n, 0, 1)
  
  x3 <- x4 <- rep(NA, m)
  
  id <- sample(1:n, m, replace = TRUE)
  
  for (g in 1:n) {
    
    x3[id == g] <- v3[g]
    x4[id == g] <- v4[g]
    
  }
  
  # transformed predictors
  u1 <- as.numeric(scale((x1 + x2)^2))
  u2 <- as.numeric(scale(cos(2*x2)))
  w3 <- as.numeric(scale(sin(2*v3)))
  w4 <- as.numeric(scale(-abs(v3 + v4)))
  
  u3 <- u4 <- rep(NA, n)
  
  for (g in 1:n) {
    
    u3[id == g] <- w3[g]
    u4[id == g] <- w4[g]
    
  }
  
  # propensity score
  if (gps_scen == "b") {
    mu_gps <- 1 + 0.5*u1 - 0.5*u2 - 0.5*u3 + 0.5*u4
  } else { # z_scen == "a"
    mu_gps <- 1 + 0.5*x1 - 0.5*x2 - 0.5*x3 + 0.5*x4
  }
  
  a <- aggregate(rnorm(n, mu_gps, sig_gps), by = list(id), mean)[,2]
  a_y <- rep(NA, n)
  
  for (g in 1:length(a))
    a_y[id == g] <- a[g]
  
  if (out_scen == "b") {
    mu_out <- -0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4 + a_y*(1 + 0.5*u1 - 0.5*u2 + 0.5*u3 - 0.5*u4)
  } else { # y_scen == "b"
    mu_out <- -0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4 + a_y*(1 + 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  }
  
  y <- rbinom(n, 1, plogis(mu_out))
  x <- cbind(x1, x2, x3, x4)
  
  # create simulation dataset
  sim <- list(y = y, a = a, x = x, id = id, a_y = a_y)
  
  return(sim)
  
}

predict_example <- function(a.vals, x, y.id, out_scen = c("a", "b")){
  
  v3 <- unique(x[,3])[order(unique(y.id))]
  v4 <- unique(x[,4])[order(unique(y.id))]
  
  # transformed predictors
  u1 <- as.numeric(scale((x[,1] + x[,2])^2))
  u2 <- as.numeric(scale(cos(2*x[,2])))
  w3 <- as.numeric(scale(sin(2*v3)))
  w4 <- as.numeric(scale(-abs(v3 + v4)))
  
  id <- unique(y.id)[order(unique(y.id))]
  u3 <- u4 <- rep(NA, nrow(x))
  
  for (g in id) {
    
    u3[y.id == g] <- w3[id == g]
    u4[y.id == g] <- w4[id == g]
    
  }
  
  u <- cbind(u1, u2, u3, u4)
  
  out <- rep(NA, length(a.vals))
  size <- unname(table(id))
  
  for(i in 1:length(a.vals)) {
    
    if (out_scen == "b") {
      mu_out <- aggregate(plogis(u %*% c(-0.75,-0.25,0.25,0.75) + 
                                   rep(a.vals[i],nrow(x))*(1 + u %*% c(0.5, -0.5, 0.5, -0.5))),
                          by = list(y.id), mean)[,2]
    } else { # y_scen == "b"
      mu_out <- aggregate(plogis(x %*% c(-0.75,-0.25,0.25,0.75) + 
                                   rep(a.vals[i],nrow(x))*(1 + x %*% c(0.5, -0.5, 0.5, -0.5))),
                          by = list(y.id), mean)[,2]
    }
    
    out[i] <- sum(size*mu_out)/sum(size)
    
  }
  
  return(out)
  
}
