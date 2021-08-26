bayes_erc <- function(s, star, y, s.id, id, w = NULL, x = NULL,
                      offset = NULL, weights = NULL, family = gaussian(),
                      a.vals = seq(min(a), max(a), length.out = 100),
                      shape = 1e-3, rate = 1e-3, scale = 1e6,
                      n.iter = 10000, n.adapt = 1000, thin = 10, 
                      h.a = 0.5, h.gamma = 0.1, span = 0.75) {
  
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
  star <- as.matrix(star)[s.id %in% id]
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
    ws <- cbind(rep(1, length(s.id)), star = star)
  } else {
    ws <- cbind(model.matrix(~ ., data.frame(w)), star = star)
  }
  
  sword <- order(s.id)
  stab <- table(s.id)
  s <- s[sword]
  star <- star[sword]
  ws <- ws[sword,]
  s.id <- s.id[sword]
  
  # when s is observed
  ws.tmp <- ws[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  
  # initialize exposures
  s.hat <- predict(lm(s.tmp ~ 0 + ., data = data.frame(ws.tmp)), newdata = data.frame(ws))
  a <- aggregate(s.hat, by = list(s.id), mean)[,2]
  a.s <- rep(a, stab)
  xa <- cbind(x, a - 8, (a - 8)^2, (a - 8)*x[,2]) # needs to be more general
  
  # data dimensions
  p <- ncol(x)
  q <- ncol(ws)
  o <- ncol(xa)
  
  # initialize parameters
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  
  # outcome stuff
  y_ <- family$linkinv(family$linkfun(y) - offset)
  gamma <- matrix(NA, nrow = n.adapt + n.iter, ncol = o)
  gamma[1,] <- coef(glm(y ~ 0 + xa, family = family, offset = offset))
  
  # the good stuff
  a.mat <- psi <- matrix(NA, nrow = floor(n.iter/thin), ncol = n)
  mhat.out <- est.mat <- var.mat <- matrix(NA, nrow = floor(n.iter/thin), ncol = length(a.vals))
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    # print(i)
    
    # sample S
    
    sig.s <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    mu.s <- (sig.s^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, mu.s, sig.s)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    a_ <- rnorm(n, a, h.a)
    xa_ <- cbind(x, a_ - 8, (a_ - 8)^2, (a_ - 8)*x[,2]) # needs to be more general
    
    log.eps <- dpois(y, family$linkinv(c(xa_%*%gamma[i - 1,]) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      dpois(y, family$linkinv(c(xa%*%gamma[i - 1,]) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
      dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)
    
    test <- log(runif(n))
    a <- ifelse(((test <= log.eps) & !is.na(log.eps)), a_, a)
    a.s <- rep(a, stab)
    
    # Sample pred parameters
    
    alpha_var <- solve(t(ws.tmp) %*% ws.tmp + diag(tau2[i - 1]/scale, q, q))
    alpha[i,] <- rmvnorm(1, alpha_var %*% t(ws.tmp) %*% s.tmp, tau2[i - 1]*alpha_var)
    
    tau2[i] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                          sum(c(s.tmp - c(ws.tmp %*% alpha[i,]))^2)/2)
    
    # Sample agg parameters
    
    omega2[i] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((s.hat - a.s)^2)/2)
    
    # Sample GPS parameters
    
    beta_var <- solve(t(x) %*% x + diag(sigma2[i - 1]/scale, p, p))
    beta[i,] <- rmvnorm(1, beta_var %*% t(x) %*% a, sigma2[i - 1]*beta_var)
    
    sigma2[i] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - c(x %*% beta[i,]))^2)/2)
    
    # Sample outcome model while cutting feedback
    
    xa <- cbind(x, a - 8, (a - 8)^2, (a - 8)*x[,2]) # needs to be more general
    gamma_ <- gamma0 <- gamma[i - 1,]
    
    for (j in 1:o) {
      
      gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma))
      
      log.eps <- sum(dfun(y, family$linkinv(c(xa %*% gamma_) + offset), log = TRUE)) -
        sum(dfun(y, family$linkinv(c(xa %*% gamma0) + offset), log = TRUE)) +
        dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma0[j], 0 , scale, log = TRUE)
      
      if ((log(runif(1)) <= log.eps) & !is.na(log.eps))
        gamma0[j] <- gamma_[j]
      else
        gamma_[j] <- gamma0[j]
      
    }
    
    gamma[i,] <- gamma_
    
    # update important things
    if (i > n.adapt & (i - n.adapt)%%thin == 0) {
      
      j <- (i - n.adapt)/thin
      a.mat[j,] <- a 
      
      xa.new.list <- lapply(a.vals, function(a.tmp, ...) {
        
        cbind(x, a.tmp - 8, (a.tmp - 8)^2, (a.tmp - 8)*x[,2]) # needs to be more general
        
      })
      
      xa.new <- rbind(xa, do.call(rbind, xa.new.list))
      x.new <- xa.new[-(1:n),1:ncol(x)]
      a.new <- rep(a.vals, each = n)
      colnames(x.new) <- colnames(x)
      colnames(xa.new) <- colnames(xa)
      
      # outcome models
      muhat.vals <- family$linkinv(c(xa.new %*% gamma[i,]))
      muhat <- muhat.vals[1:n]
      muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
      mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
      
      # exposure models
      pimod.vals <- c(x.new %*% beta[i,])
      pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2[i]))
      pihat.mat <- matrix(pihat.vals, nrow = n, ncol = length(a.vals))
      
      # pseudo-outcome
      psi[j,] <- (y_ - muhat)*(phat/pihat) + mhat
      
      # integrate
      mhat.out[j,] <- colMeans(muhat.mat)
      int.mat <- muhat.mat - matrix(rep(mhat, n), byrow = T, nrow = n)
      
      dr_out <- sapply(a.vals, dr_est, psi = psi[j,], a = a, family = gaussian(), 
                       span = span, int.mat = int.mat, se.fit = TRUE)
      
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
  estimate <- colMeans(est.mat)
  variance <- colMeans(var.mat) + (1 + 1/nrow(a.mat))*apply(est.mat, 2, var)
  
  rslt <- list(estimate = estimate, variance = variance, 
               accept.a = accept.a, accept.gamma = accept.gamma,
               mcmc = list(a.mat = a.mat, gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}
