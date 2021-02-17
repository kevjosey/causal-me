gen_data <- function(l, m, n, sig_epe = sqrt(2), sig_gps = 1,
                     gps_scen = c("a", "b"), out_scen = c("a", "b")) {
  
  if (m > n | l > n | m > l)
    stop("you stop that, you!")
  

  
  y.id <- sample(1:l, n, replace = TRUE)
  s.id <- sample(1:m, l, replace = TRUE)
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  w1 <- stats::rnorm(l, 0, 1)
  w2 <- stats::rnorm(l, 0, 1)
  w3 <- stats::rnorm(l, 0, 1)
  w4 <- stats::rnorm(l, 0, 1)
  u1 <- stats::rnorm(m, 0, 1)
  u2 <- stats::rnorm(m, 0, 1)
  
  x3 <- x4 <- rep(NA, n)
  
  for (g in 1:m) {
    
    x3[s.id == g] <- w3[g]
    x4[s.id == g] <- w4[g]
    
  }
  
  mu_gps <- 1 + 0.5*v[,1] - 0.5*v[,2] - 0.5*v[,3] + 0.5*v[,4]
  a <- rnorm(m, mu_gps, sig_gps)
  a_y <- rep(NA, n)
  a_s <- rep(NA, l)
  
  for (g in 1:m) {
    
    a_y[y.id == g] <- a[g]
    a_s[s.id == g] <- a[g]
    
  }
  
  s <- rnorm(l, a_s, sig_epe)
  
  if (out_scen == "b") {
    mu_out <- -0.75*u1 - 0.25*u2 + 0.25*u3 + 0.75*u4 + a_y*(1 + 0.5*u1 - 0.5*u2 + 0.5*u3 - 0.5*u4)
  } else { # y_scen == "b"
    mu_out <- -0.75*x1 - 0.25*x2 + 0.25*x3 + 0.75*x4 + a_y*(1 + 0.5*x1 - 0.5*x2 + 0.5*x3 - 0.5*x4)
  }
  
  y <- rbinom(n, 1, plogis(mu_out))
  
  # create simulation dataset
  sim <- list(s = s, y = y, x = x, s.id = s.id, y.id = y.id, a = a, a_y = a_y)
  
  return(sim)
  
}