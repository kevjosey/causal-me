
gps_erc <- function(s, star, y, s.id, id, family = gaussian(),
                      offset = rep(0, length(id)), w = NULL, x = NULL,
                      a.vals = seq(min(a), max(a), length.out = 100),
                      shape = 1e-3, rate = 1e-3, scale = 1e6,
                      thin = 10, n.iter = 10000, n.adapt = 1000, n.boot = 100,
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
  star <- as.matrix(star)[s.id %in% id]
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
  star <- star[sword]
  ws <- ws[sword,]
  s.id <- s.id[sword]
  
  ws.tmp <- ws[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  stab <- table(s.id)
  
  # initialize exposures
  
  s.hat <- predict(lm(s.tmp ~ 0 + ., data = data.frame(ws.tmp)), newdata = data.frame(ws))
  a <- aggregate(s.hat, by = list(s.id), mean)[,2]
  a.s <- rep(a, stab)
  # path <- seq(0, 1, length.out = 11)
  
  # dimensions
  p <- ncol(x)
  q <- ncol(ws)
  l <- sum(!is.na(s))
  
  # initialize parameters
  gamma <- array(NA, dim = c(n.adapt + n.iter, n.boot, ncol = deg.num + 1))
  # gamma <- matrix(NA, nrow = n.adapt + n.iter, ncol = o)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  a.mat <- ipw.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  a.mat[1,] <- a
  
  y_ <- family$linkinv(family$linkfun(y) - offset)
  nsa <- poly(a, degree = deg.num, raw = TRUE)
  xa <- cbind(1, nsa)
  o <- ncol(xa)
  # gamma[1,] <- coef(glm(y ~ 0 + xa, family = family, offset = offset))
  gamma[1,,] <- matrix(rep(coef(glm(y ~ 0 + xa, family = family, offset = offset)), n.boot),
                       byrow = TRUE, nrow = n.boot, ncol = o)
  ipw.mat[1,] <- ipw(a = a, x = x, beta = beta[1,], sigma2 = sigma2[1], a.vals = a.vals)
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    print(i)
    
    # sample S
    
    sig.s <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    mu.s <- (sig.s^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, mu.s, sig.s)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    
    # z.hat <- aggregate(s.hat, by = list(s.id), sum)[,2]
    # sig.a <- sqrt((1/sigma2[i - 1] + stab/omega2[i - 1])^(-1))
    # mu.a <- (sig.a^2)*(c(x %*% beta[i - 1,])/sigma2[i - 1] + z.hat/omega2[i - 1])
    # a <- a.mat[i,] <- rnorm(n, mu.a, sig.a)
    # a.s <- rep(a, stab)
    
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    a_ <- rnorm(n, a, h.a)
    xa_ <- cbind(1, predict(nsa, a_))
    ipw_ <- ipw(a = a_, x = x, beta = beta[i - 1,], sigma2 = sigma2[i - 1], a.vals = a.vals)

    log.eps <- ipw_*dfun(y, family$linkinv(c(xa_ %*% colMeans(gamma[i - 1,,])) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      ipw.mat[i - 1,]*dfun(y, family$linkinv(c(xa %*% colMeans(gamma[i - 1,,])) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
      dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)

    test <- log(runif(n))
    a <- a.mat[i,] <- ifelse(((test <= log.eps) & !is.na(log.eps)), a_, a)
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
    
    # Metropolis-Hastings
    # xa <- cbind(1, predict(nsa, a))
    # ipw.mat[i,] <- ipw(a = a, x = x, beta = beta[i,], sigma2 = sigma2[i], a.vals = a.vals)
    # gamma_ <- gamma0 <- gamma[i-1,]
    # 
    # for (j in 1:o) {
    # 
    #   gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma))
    # 
    #   log.eps <- sum(ipwmat[i,]*dfun(y, family$linkinv(c(xa %*% gamma_) + offset), log = TRUE)) -
    #     sum(ipwmat[i,]*dfun(y, family$linkinv(c(xa %*% gamma0) + offset), log = TRUE)) +
    #     dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma0[j], 0 , scale, log = TRUE)
    # 
    #   if ((log(runif(1)) <= log.eps) & !is.na(log.eps))
    #     gamma0[j] <- gamma_[j]
    #   else
    #     gamma_[j] <- gamma0[j]
    # 
    # }
    # 
    # gamma[i,] <- gamma_
    
    # uncertainty design
    xa <- cbind(1, predict(nsa, a))
    ipw.mat[i,] <- ipw(a = a, x = x, beta = beta[i,], sigma2 = sigma2[i], a.vals = a.vals)
    gmod <- rstanarm::stan_glm(y ~ 0 + ., data = data.frame(y = y, xa), family = poisson, QR = TRUE,
                               offset = offset, weights = ipw.mat[i,], iter = 2*n.boot, chains = 1, refresh = 0)
    gamma[i,,] <- do.call(cbind, split.along.dim(as.array(gmod), 3))
    
  }
  
  accept.a <- mean(apply(a.mat[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) ))
  # accept.gamma <- mean(apply(gamma[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) ))
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  gamma <- gamma[keep,,]
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  a.mat <- a.mat[keep,]
  ipw.mat <- ipw.mat[keep,]
  
  out <- lapply(1:nrow(a.mat), function(i, ...){
    
    a <- a.mat[i,]
    # n <- length(a)
    # x.new <- x[rep(1:n, length(a.vals) + 1), ]
    # a.new <- c(a, rep(a.vals, each = n))
    # colnames(x.new) <- colnames(x)
    
    # exposure model
    # pimod.vals <- c(x.new %*% beta[i,])
    # pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2[i]))
    # pihat <- pihat.vals[1:n]
    # pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
    # phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
    # phat[which(phat < 0)] <- 1e-6
    
    # outcome model
    # xa.new.list <- lapply(1:length(a.vals), function(j, ...) {
    # 
    #   cbind(1, matrix(rep(c(predict(nsa, a.vals[j])), n), byrow = TRUE, nrow = n), 
    #         dnorm(a.vals[j], x %*% beta[i,], sqrt(sigma2[i])))
    # 
    # })
    # 
    # xa.new <- do.call(rbind, xa.new.list)
    # colnames(xa.new) <- colnames(xa)
    
    # outcome models
    # muhat.vals <- family$linkinv(c(xa.new %*% gamma[i,]))
    # muhat <- muhat.vals[1:n]
    # muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
    # mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
    
    # integrate
    # phat.mat <- matrix(rep(colMeans(pihat.mat), n), byrow = T, nrow = n)
    # mhat.mat <- matrix(rep(colMeans(muhat.mat), n), byrow = T, nrow = n)
    # intfn <- (muhat.mat - mhat.mat) * phat.mat
    # int <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
    #                (intfn[,-1] + intfn[,-length(a.vals)]) / 2, 1, sum)
    
    psi <- y_*ipw.mat[i,]
    
    # dr_mod <- mgcv::gam(psi ~ s(a), family = gaussian(), data = data.frame(psi = psi, a = a))
    # dr_out <- predict(dr_mod, newdata = data.frame(a = a.vals))
    
    dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = NULL, span = span, se.fit = FALSE)

    # uncertainty design
    # dr_out <- sapply(1:n.boot, function(j, ...){
    # 
    #   muhat <- family$linkinv(c(xa.new %*% gamma[i,j,]))
    #   muhat.mat <- matrix(muhat, nrow = n, ncol = length(a.vals))
    #   mhat <- colMeans(muhat.mat)
    #   return(mhat)
    # 
    # })
    # 
    # estimate <- rowMeans(dr_out)
    # variance <- apply(dr_out, 1, var)
    # return(list(estimate = estimate, variance = variance))
    
    return(dr_out)
    
  })
  
  a.mat <- a.mat[,order(shield)]
  est.mat <- do.call(rbind, out)
  # var.mat <- do.call(rbind, lapply(out, function(arg, ...) arg$variance))
  estimate <- colMeans(est.mat)
  variance <- apply(est.mat, 2, var)
  hpdi <- apply(est.mat, 2, hpd)
  
  rslt <- list(estimate = estimate, variance = variance, hpdi = hpdi, accept.a = accept.a,
               mcmc = list(a.mat = a.mat, gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}
