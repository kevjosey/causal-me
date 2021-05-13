
bayes_erc <- function(s, star, y, s.id, id, family = gaussian(),
                      offset = rep(0, length(id)), w = NULL, x = NULL,
                      a.vals = seq(min(a), max(a), length.out = 100),
                      shape = 1e-3, rate = 1e-3, scale = 1e6,
                      thin = 10, n.iter = 10000, n.adapt = 1000,
                      h.a = 0.5, h.gamma = 0.1, deg.num = 2, span = 0.75) {
  
  dfun <- dpois
  
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
  star <- as.matrix(star)[s.id %in% id,,drop = FALSE]
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
    ws <- cbind(rep(1, length(s.id)), star = star)
  } else {
    ws <- cbind(model.matrix(~ ., data.frame(w)), star = star)
  }
  
  sword <- order(s.id)
  s <- s[sword]
  star <- star[sword,,drop = FALSE]
  ws <- ws[sword,]
  s.id <- s.id[sword]
  
  stab <- table(s.id)
  
  ws.tmp <- ws[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  
  # initialize exposures
  
  s.hat <- predict(lm(s.tmp ~ 0 + ., data = data.frame(ws.tmp)), newdata = data.frame(ws))
  a <- aggregate(s.hat, by = list(s.id), mean)[,2]
  a.s <- rep(a, stab)
  
  nsa <- ns(a, df = deg.num)
  xa <- cbind(nsa)
  
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
  amat <- ipwmat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  # path <- seq(0,1,length.out = 11)
  
  gamma[1,] <- coef(glm(y ~ 0 + xa, family = family, offset = offset))
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  ipwmat[1,] <- ipw(a = a, x = x, beta = beta[1,], sigma2 = sigma2[1], a.vals = a.vals) 
  amat[1,] <- a
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    # sample S
    
    sig <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    hat <- (sig^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, hat, sig)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    a_ <- rnorm(n, a, h.a)
    xa_ <- cbind(predict(nsa, a_))
    
    log.eps <- ipwmat[i - 1,]*dfun(y, family$linkinv(c(xa_%*%gamma[i - 1,]) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      ipwmat[i - 1,]*dfun(y, family$linkinv(c(xa%*%gamma[i - 1,]) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
      dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)

    test <- log(runif(n))
    a <- amat[i,] <- ifelse(((test <= log.eps) & !is.na(log.eps)), a_, a)
    
    xa <- cbind(predict(nsa, a))
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
    
    gamma_ <- gamma0 <- gamma[i,] <- gamma[i-1,]
    ipwmat[i,] <-  ipw(a = a, x = x, beta = beta[i,], sigma2 = sigma2[i], a.vals = a.vals) 
    
    for (j in 1:o) {

      gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma))

      log.eps <- sum(ipwmat[i,]*dfun(y, family$linkinv(c(xa %*% gamma_) + offset), log = TRUE)) -
        sum(ipwmat[i,]*dfun(y, family$linkinv(c(xa %*% gamma0) + offset), log = TRUE)) +
        dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma0[j], 0 , scale, log = TRUE)

      if ((log(runif(1)) <= log.eps) & !is.na(log.eps))
        gamma0[j] <- gamma_[j]
      else
        gamma_[j] <- gamma0[j]

    }
    
    # for (l in 1:length(path)){
    #   
    #   a_ <- path[l]*amat[i,] + (1 - path[l])*amat[i - 1,]
    #   ipw_ <- path[l]*ipwmat[i,] + (1 - path[l])*ipwmat[i - 1,]
    #   xa_ <-cbind(predict(nsa, a_))
    # 
    #   for (j in 1:o) {
    #     
    #     gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma))
    #     
    #     log.eps <- sum(ipw_*dfun(y, family$linkinv(c(xa_ %*% gamma_) + offset), log = TRUE)) -
    #       sum(ipw_*dfun(y, family$linkinv(c(xa_ %*% gamma0) + offset), log = TRUE)) +
    #       dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma0[j], 0 , scale, log = TRUE)
    #     
    #     if ((log(runif(1)) <= log.eps) & !is.na(log.eps))
    #       gamma0[j] <- gamma_[j]
    #     else
    #       gamma_[j] <- gamma0[j]
    #     
    #   }
    #   
    # }

    gamma[i,] <- gamma_
    
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
  
  # out <- lapply(1:nrow(amat), function(i, ...){
  # 
  #   a <- amat[i,]
  #   xa <- cbind(x, predict(nsa, a))
  # 
  #   xa.new.list <- lapply(a.vals, function(a.tmp, ...) {
  # 
  #     cbind(x, matrix(rep(c(predict(nsa, a.tmp)), n), byrow = TRUE, nrow = n))
  # 
  #   })
  # 
  #   xa.new <- rbind(xa, do.call(rbind, xa.new.list))
  #   x.new <- xa.new[,1:ncol(x)]
  #   a.new <- c(a, rep(a.vals, each = n))
  #   colnames(x.new) <- colnames(x)
  #   colnames(xa.new) <- colnames(xa)
  # 
  #   # exposure models
  #   pimod.vals <- c(x.new %*% beta[i,])
  #   pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2[i]))
  #   pihat <- pihat.vals[1:n]
  #   pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  #   phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
  #   phat[which(phat < 0)] <- 1e-6
  # 
  #   # outcome models
  #   muhat.vals <- family$linkinv(c(xa.new %*% gamma[i,]))
  #   muhat <- muhat.vals[1:n]
  #   muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  #   mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
  # 
  #   # integrate
  #   # phat.mat <- matrix(rep(colMeans(pihat.mat), n), byrow = T, nrow = n)
  #   # mhat.mat <- matrix(rep(colMeans(muhat.mat), n), byrow = T, nrow = n)
  #   # intfn <- (muhat.mat - mhat.mat) * phat.mat
  #   # int <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
  #   #                (intfn[,-1] + intfn[,-length(a.vals)]) / 2, 1, sum)
  # 
  #   psi <- (y.new - muhat)/(pihat/phat) + mhat
  # 
  #   # dr_mod <- mgcv::gam(psi ~ s(a), family = gaussian(), data = data.frame(psi = psi, a = a))
  #   # dr_out <- predict(dr_mod, newdata = data.frame(a = a.vals))
  # 
  #   dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = NULL, span = span, se.fit = FALSE)
  # 
  #   return(dr_out)
  # 
  # })
  
  xa.new <- cbind(predict(nsa, a.vals))
  est.mat <- sapply(1:nrow(amat), function(i, ...) {
    
    family$linkinv(xa.new%*%gamma[i,])
    
  })
  
  amat <- amat[,order(shield)]
  estimate <- rowMeans(est.mat)
  variance <- apply(est.mat, 1, var)
  hpdi <- apply(est.mat, 1, hpd)
  
  rslt <- list(estimate = estimate, variance = variance, hpdi = hpdi,
               mcmc = list(amat = amat, gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2),
               accept.a = accept.a, accept.gamma = accept.gamma)
  
  return(rslt)
  
}

ipw <- function(a, x, beta, sigma2, a.vals) {
  
  n <- length(a)
  x.new <- rbind(x, x[rep(1:n, length(a.vals)), ])
  a.new <- c(a, rep(a.vals, each = n))
  pimod.vals <- c(x.new %*% beta)
  pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2))
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
  phat[which(phat < 0)] <- 1e-6
  out <- phat/pihat
  return(out)
  
}

# highest posterior density
hpd <- function(x, alpha = 0.05){
  
  n <- length(x)
  m <- round(n * alpha)
  x <- sort(x)
  y <- x[(n - m + 1):n] - x[1:m]
  z <- min(y)
  k <- which(y == z)[1]
  c(x[k], x[n - m + k])
  
}
