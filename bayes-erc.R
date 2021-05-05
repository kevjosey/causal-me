
bayes_erc <- function(s, star, y, s.id, id, family = gaussian(),
                      offset = rep(0, length(id)), w = NULL, x = NULL,
                      a.vals = seq(min(a), max(a), length.out = 100),
                      shape = 1e-3, rate = 1e-3, scale = 1e6,
                      thin = 10, n.iter = 10000, n.adapt = 1000,
                      h.a = 0.5, h.gamma = 0.1, deg.num = 2, span = 0.75) {
  
  # remove any s.id not present in id
  check <- unique(s.id)[order(unique(s.id))]
  check <- check[check %in% id]
  
  if(any(duplicated(id)))
    stop("id must be unique.")
  
  if(length(check) < length(id))
    stop("some observations in id are not represented by measurements of s.id.")
  
  if(length(check) > length(id))
    warning("deleting some exposures without an associated outcome.")
  
  s <- s[s.id %in% id]
  star <- star[s.id %in% id]
  w <- w[s.id %in% id,]
  s.id <- s.id[s.id %in% id]
  
  # create variables
  m <- length(s.id)
  n <- length(id)
  
  if (is.null(x)) {
    x <- matrix(1, nrow = length(id), ncol = 1)
  } else {
    x <- model.matrix(~ ., data.frame(x))
  }
  
  shield <- order(id)
  y <- y[shield]
  x <- x[shield,]
  id <- id[shield]
  
  if (is.null(w)) {
    # ws <- cbind(rep(1, length(s.id)), star)
    ws <- cbind(rep(1, length(s.id)), poly(star, degree = deg.num))
  } else {
    # ws <- cbind(model.matrix(~ ., data.frame(w)), star = star)
    ws <- cbind(model.matrix(~ ., data.frame(w)), star = poly(star, degree = deg.num))
  }
  
  sword <- order(s.id)
  s <- s[sword]
  star <- star[sword]
  ws <- ws[sword,]
  s.id <- s.id[sword]
  
  stab <- table(s.id)
  
  ws.tmp <- ws[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  
  # initialize exposures
  
  s.hat <- predict(lm(s.tmp ~ 0 + ., data = data.frame(ws.tmp)), newdata = data.frame(ws))
  a <- aggregate(s.hat, by = list(s.id), mean)[,2]
  a.s <- rep(a, stab)
  
  nsa <- poly(a, degree = deg.num)
  xa <- as.matrix(cbind(x, nsa))
  
  # dimensions
  p <- ncol(x)
  q <- ncol(ws)
  o <- ncol(xa)
  l <- sum(!is.na(s))
  
  # initialize parameters
  gamma <- matrix(NA, nrow = n.adapt + n.iter, ncol = o)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  
  gamma[1,] <- coef(glm(y ~ 0 + xa, family = family, offset = offset))
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(star - a.s)
  
  amat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    # sample S
    
    sig <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    hat <- (sig^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, hat, sig)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    
    a_ <-  rnorm(n, a, h.a)
    xa_ <- as.matrix(cbind(x, predict(nsa, a_)))
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    
    test <- log(runif(n))
    
    log.eps <- dpois(y, exp(c(xa_%*%gamma[i - 1,]) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) + 
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      dpois(y, exp(c(xa%*%gamma[i - 1,]) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
      dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)
    
    a <- amat[i,] <- ifelse(((test <= log.eps) & !is.na(log.eps)), a_, a)
    
    xa <- as.matrix(cbind(x, predict(nsa, a)))
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
    
    gamma_ <- gamma[i,] <- gamma[i-1,]
    
    for (j in 1:o){
      
      gamma_[j] <- c(rnorm(1, gamma[i,j], h.gamma))
      
      log.eps <- sum(dpois(y, exp(c(xa %*% gamma_) + offset), log = TRUE)) -
        sum(dpois(y, exp(c(xa %*% gamma[i,]) + offset), log = TRUE)) +
        dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma[i,j], 0 , scale, log = TRUE)
      
      if ((log(runif(1)) <= log.eps) & !is.na(log.eps))
        gamma[i,j] <- gamma_[j]
      else
        gamma_[j] <- gamma[i,j]
      
    }
    
  }
  
  accept.a <- mean(apply(amat[(n.adapt + 1):nrow(amat),], 2, function(x) mean(diff(x) != 0) ))
  accept.gamma <- mean(apply(gamma[(n.adapt + 1):nrow(gamma),], 2, function(x) mean(diff(x) != 0) ))
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  gamma <- gamma[keep,]
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  amat <- amat[keep,]
  
  y.new <- exp(log(y) - offset)
  
  out <- lapply(1:nrow(amat), function(i, ...){
    
    a <- amat[i,]
    xa <- as.matrix(cbind(x, predict(nsa, a)))
    
    if (ncol(x) > 1) {
      pimod.vals <- c(x %*% beta[i,])
    } else {
      pimod.vals <- c(x %*% beta[i])
    }
    
    # exposure models
    pihat <- dnorm(a, pimod.vals, sqrt(sigma2[i]))
    phat.vals <- sapply(a.vals, function(a.tmp, ...)
      mean(dnorm(a.tmp, pimod.vals, sqrt(sigma2[i]))))
    phat <- predict(smooth.spline(a.vals, phat.vals), x = a)$y
    phat[which(phat < 0)] <- 1e-6
    # phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)
    
    muhat <- exp(c(xa %*% gamma[i,]))
    
    # predict marginal outcomes given a.vals (or a.agg)
    muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
      
      xa.tmp <- as.matrix(cbind(x, matrix(rep(c(predict(nsa, a.tmp)), n), byrow = TRUE, nrow = n)))
      return(exp(c(xa.tmp %*% gamma[i,])))
      
    })
    
    # aggregate muhat.vals and integrate for influence function
    mhat.vals <- colMeans(muhat.mat)
    mhat <- predict(smooth.spline(a.vals, mhat.vals), x = a)$y
    # mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
    
    # integrate
    # intfn <- (muhat.mat - mhat.mat) * phat.mat
    # int <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
    #                (intfn[,-1] + intfn[,-length(a.vals)]) / 2, 1, sum)
    
    psi <- c((y.new - muhat)/(pihat/phat) + mhat)
    
    # dr_mod <- loess(psi ~ a, data = data.frame(psi = psi, a = a), span = span, degree = 1)
    # dr_out <- predict(dr_mod, newdata = data.frame(a = a.vals))
    
    dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = NULL,
                     family = family, span = span, se.fit = FALSE)
    
    return(dr_out)
    
  })
  
  amat <- amat[,order(shield)]
  est.mat <- do.call(rbind, out)
  estimate <- colMeans(est.mat)
  variance <- apply(est.mat, 2, var)
  
  rslt <- list(estimate = estimate, variance = variance,
               mcmc = list(gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2, amat = amat),
               accept.a = accept.a, accept.gamma = accept.gamma)
  
  return(rslt)
  
}
