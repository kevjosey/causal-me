# WARNING: Only implements DR Estimation which violates congeniality concepts
bayes_erf <- function(s, t, y, s.id, id, w = NULL, x = NULL, df = 4,
                      offset = NULL, weights = NULL, family = gaussian(),
                      a.vals = seq(min(s), max(s), length.out = 100),
                      shape = 1e-3, rate = 1e-3, scale = 1e6,
                      n.iter = 10000, n.adapt = 1000, thin = 10, 
                      h.a = 0.5, span = 0.75) {
  
  # remove any s.id not present in id
  check <- unique(s.id)[order(unique(s.id))]
  check <- check[check %in% id]
  
  if(any(duplicated(id)))
    stop("id must be unique.")
  
  if(length(check) < length(id))
    stop("some observations in id are not represented by measurements of s.id.")
  
  if(length(check) > length(id))
    warning("deleting some exposures without an associated outcome.")
  
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  s <- s[s.id %in% id]
  t <- as.matrix(t)[s.id %in% id]
  w <- w[s.id %in% id,]
  s.id <- s.id[s.id %in% id]
  
  # create variables
  m <- length(s.id)
  n <- length(id)
  l <- sum(!is.na(s))
  
  if (is.null(x)) {
    x <- matrix(1, nrow = length(id), ncol = 1)
  } else {
    x <- model.matrix(~ ., data.frame(x))
  }
  
  shield <- order(id)
  y <- y[shield]
  x <- x[shield,]
  id <- id[shield]
  weights <- weights[shield]
  offset <- offset[shield]
  
  if (is.null(w)) {
    ws <- cbind(rep(1, length(s.id)), t = t)
  } else {
    ws <- cbind(model.matrix(~ ., data.frame(w)), t = t)
  }
  
  sword <- order(s.id)
  stab <- table(s.id)
  s <- s[sword]
  t <- t[sword]
  ws <- ws[sword,]
  s.id <- s.id[sword]
  
  # when s is observed
  ws.tmp <- ws[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  
  # initialize exposures
  s.hat <- predict(lm(s.tmp ~ 0 + ., data = data.frame(ws.tmp)), newdata = data.frame(ws))
  a <- aggregate(s.hat, by = list(s.id), mean)[,2]
  a.s <- rep(a, stab)
  nsa <- ns(a, df = df)
  xa <- cbind(x, nsa, a*x[,2]) # needs to be more general
  
  # data dimensions
  p <- ncol(x)
  q <- ncol(ws)
  o <- ncol(xa)
  
  # initialize parameters
  gamma <- matrix(NA, nrow = n.adapt + n.iter, ncol = o)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  
  fit.a <- lm(a ~ 0 + x)
  fit.s <- lm(s.tmp ~ 0 + ws.tmp)
  fit.y <- glm(y ~ 0 + xa, family = family, offset = offset)
  
  gamma[1,] <- coef(fit.y)
  beta[1,] <- coef(fit.a)
  alpha[1,] <- coef(fit.s)
  sigma2[1] <- sigma(fit.a)^2
  tau2[1] <- sigma(fit.s)^2
  omega2[1] <- var(s.hat - a.s)
  
  # auxilliary
  ybar <- family$linkinv(family$linkfun(y) - offset)
  h.gamma <- 0.5*sqrt(diag(vcov(fit.y)))
  h.a <- rep(h.a, n)
  accept.a <- rep(0, n)
  accept.gamma <- rep(0, o)
  
  # the good stuff
  a.mat <- psi <- matrix(NA, nrow = floor(n.iter/thin), ncol = n)
  mhat.out <- est.mat <- var.mat <- matrix(NA, nrow = floor(n.iter/thin), ncol = length(a.vals))
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    print(i)
    
    # sample S
    sig.s <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    mu.s <- (sig.s^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, mu.s, sig.s)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    a_ <- rnorm(n, a, h.a)
    xa_ <- cbind(x, predict(nsa, newx = a_), a_*x[,2]) # needs to be more general
    
    log.eps <- dpois(y, family$linkinv(c(xa_%*%gamma[i - 1,]) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      dpois(y, family$linkinv(c(xa%*%gamma[i - 1,]) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
      dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)
    
    test <- log(runif(n)) <= log.eps & !is.na(log.eps)
    a <- ifelse(test, a_, a)
    accept.a <- accept.a + test
    a.s <- rep(a, stab)
    
    # sample pred parameters
    alpha_var <- solve(t(ws.tmp) %*% ws.tmp + diag(tau2[i - 1]/scale, q, q))
    alpha[i,] <- rmvnorm(1, alpha_var %*% t(ws.tmp) %*% s.tmp, tau2[i - 1]*alpha_var)
    tau2[i] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                          sum(c(s.tmp - c(ws.tmp %*% alpha[i,]))^2)/2)
    
    # sample agg parameters
    omega2[i] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((s.hat - a.s)^2)/2)
    
    # sample gps parameters
    beta_var <- solve(t(x) %*% x + diag(sigma2[i - 1]/scale, p, p))
    beta[i,] <- rmvnorm(1, beta_var %*% t(x) %*% a, sigma2[i - 1]*beta_var)
    sigma2[i] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - c(x %*% beta[i,]))^2)/2)
    
    # sample outcome model while cutting feedback
    xa <- cbind(x, predict(nsa, newx = a), a*x[,2]) # needs to be more general
    gamma_ <- gamma0 <- gamma[i - 1,]
    
    for (j in 1:o) {
      
      gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma[j]))
      
      log.eps <- sum(dpois(y, family$linkinv(c(xa %*% gamma_) + offset), log = TRUE)) -
        sum(dpois(y, family$linkinv(c(xa %*% gamma0) + offset), log = TRUE)) +
        dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma0[j], 0, scale, log = TRUE)
      
      exam <- (log(runif(1)) <= log.eps) & !is.na(log.eps)
      accept.gamma[j] <- accept.gamma[j] + exam
      
      if (exam)
        gamma0[j] <- gamma_[j]
      else
        gamma_[j] <- gamma0[j]
      
    }
    
    gamma[i,] <- gamma_
    
    if(ceiling(i/100)==floor(i/100) & i < n.adapt) {
    
      h.gamma <- ifelse(accept.gamma > 40, h.gamma + 0.1 * h.gamma,
                        ifelse(accept.gamma < 30, h.gamma - 0.1 * h.gamma, h.gamma))
      
      h.a <- ifelse(accept.a > 20, h.a + 0.1 * h.a,
                    ifelse(accept.a < 10, h.a - 0.1 * h.a, h.a))
      
      accept.a <- rep(0, n)
      accept.gamma <- rep(0, o)
      
    }
    
    # update important things
    if (i > n.adapt & (i - n.adapt)%%thin == 0) {
      
      j <- (i - n.adapt)/thin
      a.mat[j,] <- a
      nsa.vals <- predict(nsa, newx = a.vals)
      
      # outcome model
      muhat <- family$linkinv(c(xa %*% gamma[i,]))
      
      muhat.mat <- sapply(a.vals, function(a.tmp, ...){
        xa.tmp <- cbind(x, nsa.vals[rep(which(a.vals == a.tmp), nrow(x)),], a.tmp*x[,2])
        family$linkinv(c(xa.tmp %*% gamma[i,]))
      })
      
      mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
      
      # exposure model for integtion
      a.std <- c(c(a, a.vals) - mean(a)) / sd(a)
      dens <- density(a.std[1:n])
      phat.vals <- approx(x = dens$x, y = dens$y, xout = a.std)$y / sd(a)
      phat <- phat.vals[1:n]
      pimod.vals <- c(x %*% beta[i,])
      a.denom <- c(a - pimod.vals) / sqrt(sigma2[i])
      dens.denom <- density(a.denom)
      pihat <- approx(x = dens.denom$x, y = dens.denom$y, xout = a.denom)$y / sqrt(sigma2[i])
      
      # integration matrix
      phat.mat <- matrix(rep(phat.vals[-(1:n)], n), byrow = T, nrow = n)
      mhat.mat <- matrix(rep(colMeans(muhat.mat), n), byrow = T, nrow = n)
      int.mat <- (muhat.mat - mhat.mat) * phat.mat
      
      # pseudo-outcome
      psi[j,] <- c(ybar - muhat)*(phat/pihat) + mhat
      
      dr_out <- sapply(a.vals, loess_est, psi = psi[j,], a = a, span = span,
                       family = poisson(), weights = weights, int.mat = int.mat, 
                       se.fit = TRUE)
      
      est.mat[j,] <- dr_out[1,]
      var.mat[j,] <- dr_out[2,]
      
    }
    
  }
  
  accept.a <- apply(a.mat, 2, function(x) mean(diff(x) != 0) )
  accept.gamma <- apply(gamma[(n.adapt + 1):(n.iter + n.adapt),], 2, function(x) mean(diff(x) != 0) )
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  gamma <- gamma[keep,]
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  
  a.mat <- a.mat[,order(shield)]
  dr_estimate <- colMeans(est.mat)
  dr_variance <- colMeans(var.mat) + (1 + 1/nrow(a.mat))*apply(est.mat, 2, var)
  
  rslt <- list(dr_estimate = dr_estimate, dr_variance = dr_variance, 
               accept.a = accept.a, accept.gamma = accept.gamma,
               mcmc = list(a.mat = a.mat, beta = beta, alpha = alpha, gamma = gamma,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}
