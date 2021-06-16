
mi_erc <- function(s, star, y, s.id, id, family = gaussian(),
                   offset = rep(0, length(id)), w = NULL, x = NULL,
                   sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.earth"),
                   a.vals = seq(min(a), max(a), length.out = 100),deg.num = 2, span = 0.75,
                   shape = 1e-3, rate = 1e-3, scale = 1e6, h.a = 0.5, h.gamma = 0.1,
                   thin = 10, n.iter = 10000, n.adapt = 100) {
  
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
  star <- star[s.id %in% id]
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
  
  # sort on id
  shield <- order(id)
  y <- y[shield]
  x <- x[shield,]
  id <- id[shield]
  
  if (is.null(w)) {
    ws <- cbind(rep(1, length(s.id)), star = star)
  } else {
    ws <- cbind(model.matrix(~ ., data.frame(w)), star = star)
  }
  
  # sort on s.id
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
  xa <- cbind(x, a - 8, (a - 8)^2, x[,2]*(a - 8))
  
  # dimensions
  p <- ncol(x)
  q <- ncol(ws)
  o <- ncol(xa)
  
  if (length(h.gamma) == 1)
    h.gamma <- rep(h.gamma, ncol(xa))
  
  if (length(h.gamma) != ncol(xa))
    stop("length(h.gamma) != ncol(xa))")
  
  # initialize parameters
  gamma <- matrix(NA, nrow = n.adapt + n.iter, ncol = o)
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  a.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  
  gamma[1,] <- coef(glm(y ~ 0 + xa, family = family, offset = offset))
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  a.mat[1,] <- a
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
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
    xa_ <- cbind(x, a_ - 8, (a_ - 8)^2, x[,2]*(a_ - 8))
    colnames(xa_) <- colnames(xa)

    log.eps <- dpois(y, family$linkinv(c(xa_ %*% gamma[i - 1,]) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      dpois(y, family$linkinv(c(xa %*% gamma[i - 1,]) + offset), log = TRUE) -
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
    
    # sample outcome
    
    xa <- cbind(x, a - 8, (a - 8)^2, x[,2]*(a - 8))
    gamma_ <- gamma0 <- gamma[i-1,]

    for (j in 1:o) {

      gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma[j]))

      log.eps <- sum(dfun(y, family$linkinv(c(xa %*% gamma_) + offset), log = TRUE)) -
        sum(dfun(y, family$linkinv(c(xa %*% gamma0) + offset), log = TRUE)) +
        dnorm(gamma_[j], 0, scale, log = TRUE) - dnorm(gamma0[j], 0 , scale, log = TRUE)

      if ((log(runif(1)) <= log.eps) & !is.na(log.eps))
        gamma0[j] <- gamma_[j]
      else
        gamma_[j] <- gamma0[j]

    }

    gamma[i,] <- gamma_
    
  }
  
  accept.a <- apply(a.mat[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) )
  accept.gamma <- apply(gamma[(n.adapt + 1):nrow(gamma),], 2, function(x) mean(diff(x) != 0) )
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  gamma <- gamma[keep,]
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  a.mat <- a.mat[keep,]
  
  a.list <- as.list(data.frame(t(a.mat)))
  
  dr_out <- lapply(a.list, erc, y = y, x = x[,-1], offset = offset, deg.num = deg.num, 
                   a.vals = a.vals, family = family, span = span, sl.lib = sl.lib)
  
  a.mat <- a.mat[,order(shield)]
  est.mat <- do.call(rbind, lapply(dr_out, function(arg, ...) arg$estimate))
  var.mat <- do.call(rbind, lapply(dr_out, function(arg, ...) arg$variance))
  estimate <- colMeans(est.mat)
  variance <- colMeans(var.mat) + (1 + 1/nrow(a.mat))*apply(est.mat, 2, var)
  
  rslt <- list(estimate = estimate, variance = variance, 
               accept.a = accept.a, accept.gamma = accept.gamma,
               mcmc = list(a.mat = a.mat, gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}

