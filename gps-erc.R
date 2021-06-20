
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
  
  # data dimensions
  p <- ncol(x)
  q <- ncol(ws)
  o <- deg.num + 3
  
  # initialize parameters
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  a.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  a.mat[1,] <- a
  
  # outcome stuff
  y_ <- family$linkinv(family$linkfun(y) - offset)
  gps <- dnorm(a, c(x%*%beta[1,]), sqrt(sigma2[1]))
  xa <- cbind(model.matrix( ~ .^2 , data.frame(gps = gps - mean(gps), a = a - mean(a))), 
              poly(a - mean(a), degree = deg.num, raw = TRUE)[,-1])
  # gamma <- matrix(NA, nrow = n.adapt + n.iter, ncol = o)
  # gamma[1,] <- coef(glm(y ~ 0 + xa, family = family, offset = offset))
  gamma <- array(NA, dim = c(n.adapt + n.iter, n.boot, ncol = o))
  gamma[1,,] <- matrix(rep(coef(glm(y ~ 0 + xa, family = family, offset = offset)), n.boot),
                       byrow = TRUE, nrow = n.boot, ncol = o)
  
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
    gps_ <- dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]))
    xa_ <- cbind(model.matrix( ~ .^2 , data.frame(gps = gps_ - mean(gps_), a = a_ - mean(a_))), 
                poly(a_ - mean(a_), degree = deg.num, raw = TRUE)[,-1])
    
    log.eps <- dfun(y, family$linkinv(c(xa_%*%colMeans(gamma[i - 1,,])) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      dfun(y, family$linkinv(c(xa%*%colMeans(gamma[i - 1,,])) + offset), log = TRUE) -
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
    # gps <- dnorm(a, c(x%*%beta[i,]), sqrt(sigma2[i]), log = TRUE)
    # xa <- cbind(model.matrix( ~ .^2 , data.frame(gps = gps - mean(gps), a = a - mean(a))), 
    #             poly(a - mean(a), degree = deg.num, raw = TRUE)[,-1],
    #             poly(gps - mean(gps), degree = deg.num, raw = TRUE)[,-1])
    # 
    # gamma_ <- gamma0 <- gamma[i - 1,]
    # 
    # for (j in 1:o) {
    #   
    #   gamma_[j] <- c(rnorm(1, gamma0[j], h.gamma))
    #   
    #   log.eps <- sum(dfun(y, family$linkinv(c(xa %*% gamma_) + offset), log = TRUE)) -
    #     sum(dfun(y, family$linkinv(c(xa %*% gamma0) + offset), log = TRUE)) +
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
    gps <- dnorm(a, c(x%*%beta[i,]), sqrt(sigma2[i]))
    xa <- cbind(model.matrix( ~ .^2 , data.frame(gps = gps - mean(gps), a = a - mean(a))),
                poly(a - mean(a), degree = deg.num, raw = TRUE)[,-1])
    gmod <- rstanarm::stan_glm(y ~ 0 + ., data = data.frame(y = y, xa), family = poisson, QR = TRUE,
                               offset = offset, iter = 2*n.boot, chains = 1, refresh = 0)
    gamma[i,,] <- do.call(cbind, split.along.dim(as.array(gmod), 3))
    
  }
  
  accept.a <- apply(a.mat[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) )
  # accept.gamma <- apply(gamma[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) )
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  gamma <- gamma[keep,,]
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  a.mat <- a.mat[keep,]

  out <- lapply(1:nrow(a.mat), function(i, ...){

    a <- a.mat[i,]
    n <- length(a)
    x.new <- x[rep(1:n, length(a.vals)), ]
    a.new <- rep(a.vals, each = n)
    # sig.a <- sqrt((1/sigma2[i] + rep(stab, length(a.vals))/omega2[i])^(-1))
    # mu.a <- (sig.a^2)*(c(x.new %*% beta[i,])/sigma2[i] + rep(z.mat[i,], length(a.vals))/omega2[i])
    pimod.vals <- c(x.new %*% beta[i,])
    pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2[i]))
    pihat.mat <- matrix(pihat.vals, nrow = n, ncol = length(a.vals))
    phat <- colMeans(pihat.mat)

    # uncertainty design
    dr_out <- sapply(1:length(a.vals), function(j, ...){

      a.std <- a - a.vals[j]
      k <- floor(min(span, 1)*length(a))
      idx <- order(abs(a.std))[1:k]
      a.std <- a.std[idx]
      max.a.std <- max(abs(a.std))
      y_ <- y_[idx]
      k.std <- c((1 - abs(a.std/max.a.std)^3)^3)*mean(pihat.mat[idx,j])/pihat.mat[idx,j]
      gh <- cbind(1, a.std)
      b <- optim(par = c(0,0), fn = opt_fun, k.std = k.std, psi = y_, gh = gh, family = family)
      mu <- family$linkinv(c(b$par[1]))

      gh.inv <- solve(t(gh) %*% diag(k.std) %*% gh)
      v.inf <- (y_ - family$linkinv(c(gh%*%b$par)))^2
      sig <- gh.inv %*% t(gh) %*% diag(k.std) %*% diag(v.inf) %*% diag(k.std) %*% gh %*% gh.inv
      return(c(mu = mu, sig = sig[1,1]))

    })

    estimate <- dr_out[1,]
    variance <- dr_out[2,]

    return(list(estimate = estimate, variance = variance))

  })
  
  a.mat <- a.mat[,order(shield)]
  est.mat <- do.call(rbind, lapply(out, function(arg, ...) arg$estimate))
  var.mat <- do.call(rbind, lapply(out, function(arg, ...) arg$variance))
  estimate <- colMeans(est.mat)
  variance <- colMeans(var.mat) + (1 + 1/nrow(a.mat))*apply(est.mat, 2, var)
  hpdi <- apply(est.mat, 2, hpd)
  
  rslt <- list(estimate = estimate, variance = variance, hpdi = hpdi,
               mcmc = list(a.mat = a.mat, gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}
