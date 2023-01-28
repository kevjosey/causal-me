
bayes_erf <- function(s, s.tilde, y, s.id, id, w = NULL, x = NULL, offset = NULL,
                      dr = FALSE, a.vals = seq(min(s), max(s), length.out = 100),
                      shape = 1e-3, rate = 1e-3, scale = 1e6,
                      n.iter = 10000, n.adapt = 1000, thin = 10,
                      bw = NULL, bw.seq = seq(0.1, 2, by = 0.1), folds = 5) {
  
  # remove any s.id not present in id
  check <- unique(s.id)[order(unique(s.id))]
  check <- check[check %in% id]
  
  if(any(duplicated(id)))
    stop("id must be unique.")
  
  if(length(check) < length(id))
    stop("some observations in id are not represented by measurements of s.id.")
  
  if(length(check) > length(id))
    warning("deleting some exposures without an associated outcome.")
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  weights <- exp(offset)
  
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
  x <- x[shield,,drop = FALSE]
  id <- id[shield]
  offset <- offset[shield]
  weights <- weights[shield]
  
  if (is.null(w)) {
    ws <- cbind(rep(1, length(s.id)), s.tilde = s.tilde)
  } else {
    ws <- cbind(model.matrix(~ ., data.frame(w)), s.tilde = s.tilde)
  }
  
  s <- s[s.id %in% id]
  s.tilde <- s.tilde[s.id %in% id]
  ws <- ws[s.id %in% id,,drop = FALSE]
  s.id <- s.id[s.id %in% id]
  
  sword <- order(s.id)
  stab <- table(s.id)
  s <- s[sword]
  s.tilde <- s.tilde[sword]
  ws <- ws[sword,,drop = FALSE]
  s.id <- s.id[sword]
  
  # when s is observed
  ws.obs <- ws[!is.na(s),,drop = FALSE]
  s.obs <- s[!is.na(s)]
  
  # initialize exposures
  s.hat <- predict(lm(s.obs ~ 0 + ., data = data.frame(ws.obs)), newdata = data.frame(ws))
  a <- aggregate(s.hat, by = list(s.id), mean)[,2]
  a.s <- rep(a, stab)
  xa <- cbind(x, a - 10, cos(pi*(a - 6)/4), (a - 10)*x[,2]) # needs to be more general
  
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
  fit.s <- lm(s.obs ~ 0 + ws.obs)
  fit.y <- glm(y ~ 0 + xa, family = poisson(), offset = offset)
  
  gamma[1,] <- coef(fit.y)
  beta[1,] <- coef(fit.a)
  alpha[1,] <- coef(fit.s)
  sigma2[1] <- sigma(fit.a)^2
  tau2[1] <- sigma(fit.s)^2
  omega2[1] <- var(s.hat - a.s)
  
  # auxilliary
  ybar <- exp(log(y) - offset)
  h.a <- rep(0.5*sd(a), n)
  h.gamma <- 0.5*sqrt(diag(vcov(fit.y)))
  accept.a <- rep(0, n)
  accept.gamma <- rep(0, o)
  
  # the good stuff
  a.mat <- psi.mat <- matrix(NA, nrow = floor(n.iter/thin), ncol = n)
  est.mat <- var.mat <- matrix(NA, nrow = floor(n.iter/thin), ncol = length(a.vals))
  
  if (dr) {
    
    est.dr <- var.dr <- matrix(NA, nrow = floor(n.iter/thin), ncol = length(a.vals))
    
  }
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    print(i)
    
    # sample S
    sig.s <- (1/omega2[i - 1] + 1/tau2[i - 1])^(-1)
    mu.s <- (sig.s)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, mu.s, sqrt(sig.s))
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    a_ <- rnorm(n, a, h.a)
    xa_ <- cbind(x, a_ - 10, cos(pi*(a_ - 6)/4), (a_ - 10)*x[,2]) # needs to be more general
    
    log.eps <- dpois(y, exp(c(xa_%*%gamma[i - 1,]) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(z.hat, a_, sqrt(omega2[i - 1]/stab), log = TRUE) -
      dpois(y, exp(c(xa%*%gamma[i - 1,]) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
      dnorm(z.hat, a, sqrt(omega2[i - 1]/stab), log = TRUE)
    
    test <- log(runif(n)) <= log.eps & !is.na(log.eps)
    a <- ifelse(test, a_, a)
    accept.a <- accept.a + test
    a.s <- rep(a, stab)
    
    # sample pred parameters
    alpha.var <- solve(t(ws.obs) %*% ws.obs + diag(tau2[i - 1]/scale, q, q))
    alpha[i,] <- rmvnorm(1, alpha.var %*% t(ws.obs) %*% s.obs, tau2[i - 1]*alpha.var)
    tau2[i] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                          sum(c(s.obs - c(ws.obs %*% alpha[i,]))^2)/2)
    
    # sample agg parameters
    omega2[i] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((s.hat - a.s)^2)/2)
    
    # sample gps parameters
    beta.var <- solve(t(x) %*% x + diag(sigma2[i - 1]/scale, p, p))
    beta[i,] <- rmvnorm(1, beta.var %*% t(x) %*% a, sigma2[i - 1]*beta.var)
    sigma2[i] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - c(x %*% beta[i,]))^2)/2)
    
    # sample outcome model while cutting feedback
    xa <- cbind(x, a - 10, cos(pi*(a - 6)/4), (a - 10)*x[,2])
    gamma_ <- gamma0 <- gamma[i - 1,]
    
    for (j in 1:o) {
      
      gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma[j]))
      
      log.eps <- sum(dpois(y, exp(c(xa %*% gamma_) + offset), log = TRUE)) -
        sum(dpois(y, exp(c(xa %*% gamma0) + offset), log = TRUE)) +
        dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma0[j], 0, scale, log = TRUE)
      
      exam <- log(runif(1)) <= log.eps & !is.na(log.eps)
      accept.gamma[j] <- accept.gamma[j] + exam
      
      if (exam)
        gamma0[j] <- gamma_[j]
      else
        gamma_[j] <- gamma0[j]
      
    }
    
    gamma[i,] <- gamma_
    
    # update random walk parameters
    if(ceiling(i/100) == floor(i/100) & i < n.adapt) {
    
      h.gamma <- ifelse(accept.gamma > 40, h.gamma + 0.1 * h.gamma,
                        ifelse(accept.gamma < 20, h.gamma - 0.1 * h.gamma, h.gamma))
      
      h.a <- ifelse(accept.a > 40, h.a + 0.1 * h.a,
                    ifelse(accept.a < 20, h.a - 0.1 * h.a, h.a))
      
      accept.a <- rep(0, n)
      accept.gamma <- rep(0, o)
      
    }
    
    # update important things
    if (i > n.adapt & (i - n.adapt)%%thin == 0) {
      
      j <- (i - n.adapt)/thin
      a.mat[j,] <- a
      
      # outcome model
      muhat <- exp(c(xa %*% gamma[i,]))
      
      muhat.mat <- sapply(a.vals, function(a.tmp, ...){
        
        xa.tmp <- cbind(x, a.tmp - 10, cos(pi*(a.tmp - 6)/4), (a.tmp - 10)*x[,2])
        exp(c(xa.tmp %*% gamma[i,]))
        
      })
      
      mhat.vals <- colMeans(muhat.mat)
      mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), x = a)$y
      mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
      
      # pseudo outcome
      psi.mat[j,] <- psi <- c(ybar - muhat + mhat)
      
      # integration matrix
      pimod.vals <- c(x %*% beta[i,])
      
      # density estimation
      phat.vals <- sapply(a.vals, function(a.tmp, ...){
        mean(dnorm(a.tmp, pimod.vals, sqrt(sigma2[i])))
      })

      phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)  
      int.mat <- (muhat.mat - mhat.mat)*phat.mat
      
      # select bw if NULL
      if (j == 1 & is.null(bw))
        bw <- cv_bw(a = a, psi = psi[j,], folds = folds, bw.seq = bw.seq)
      
      # asymptotics
      out <- sapply(a.vals, kern_est, psi = psi, a = a, weights = weights, 
                    bw = bw, se.fit = TRUE, int.mat = int.mat, a.vals = a.vals)
      
      est.mat[j,] <- out[1,]
      var.mat[j,] <- out[2,]
      
      # double-robust model
      if (dr) {
        
        pihat <- dnorm(a, pimod.vals, sqrt(sigma2[i]))
        phat <- predict(smooth.spline(x = a.vals, y = phat.vals), x = a)$y
        phat[phat <= 0] <- .Machine$double.eps
        psi.dr <- c(ybar - muhat)*(phat/pihat) + mhat

        # asymptotics
        out.dr <- sapply(a.vals, kern_est, psi = psi.dr, a = a, weights = weights,
                         bw = bw, se.fit = TRUE, int.mat = int.mat, a.vals = a.vals)

        est.dr[j,] <- out.dr[1,]
        var.dr[j,] <- out.dr[2,]

      }
      
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
  psi.mat <- psi.mat[,order(shield)]
  estimate <- colMeans(est.mat)
  variance <- colMeans(var.mat) + (1 + 1/nrow(a.mat))*apply(est.mat, 2, var)
  
  rslt <- list(estimate = estimate, variance = variance, 
               accept.a = accept.a, accept.gamma = accept.gamma, 
               mcmc = list(a.mat = a.mat, psi.mat = psi.mat,
                           beta = beta, alpha = alpha, gamma = gamma,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  if (dr) {
    
    estimate.dr <- colMeans(est.dr)
    variance.dr <- colMeans(var.dr) + (1 + 1/nrow(a.mat))*apply(est.dr, 2, var)
    
    rslt <- list(estimate = estimate, variance = variance, 
                 estimate.dr = estimate.dr, variance.dr = variance.dr,
                 accept.a = accept.a, accept.gamma = accept.gamma,
                 mcmc = list(a.mat = a.mat, psi.mat = psi.mat,
                             beta = beta, alpha = alpha, gamma = gamma,
                             sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
    
  }
  
  return(rslt)
  
}
