
mi_erc <- function(s, star, y, s.id, id, family = gaussian(),
                   offset = rep(0, length(id)), w = NULL, x = NULL,
                   sl.lib = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.ranger", "SL.earth"),
                   a.vals = seq(min(a), max(a), length.out = 100),
                   shape = 1e-3, rate = 1e-3, scale = 1e6,
                   thin = 10, n.iter = 10000, n.adapt = 100, n.boot = 100,
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
  star <- as.matrix(star)[s.id %in% id,]
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

  # dimensions
  p <- ncol(x)
  q <- ncol(ws)
  l <- sum(!is.na(s))
  
  # initialize parameters
  gamma <- array(NA, dim = c(n.adapt + n.iter, n.boot, deg.num + 2))
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  amat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  muhat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  est.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = length(a.vals))
  var.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = length(a.vals))
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  amat[1,] <- a
  
  # gps <- dnorm(a, x%*%beta[1,], sqrt(sigma2[1]))
  nsa <- poly(a, degree = deg.num)
  xa <- cbind(x, nsa)
  df <- data.frame(y = y, xa[,-1])
  o <- ncol(xa)
  
  mumod <- glm(y ~ ., data = df, family = family, offset = offset)
  muhat[1,] <- predict(mumod, newdata = df, type = "response")
  y_ <- family$linkinv(family$linkfun(y) - offset)
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    print(i)
    
    # sample S
    
    sig <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    hat <- (sig^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, hat, sig)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    
    # z.hat <- aggregate(s.hat, by = list(s.id), sum)[,2]
    # 
    # sig <- sqrt((1/sigma2[i - 1] + stab/omega2[i - 1])^(-1))
    # hat <- (sig^2)*(c(x %*% beta[i - 1,])/sigma2[i - 1] + z.hat/omega2[i - 1])
    # a <- amat[i,] <- rnorm(n, hat, sig)
    # 
    # xa <- as.matrix(cbind(x, predict(nsa, a)))
    # a.s <- rep(a, stab)
    
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    a_ <- rnorm(n, a, h.a)
    xa_ <- data.frame(x, predict(nsa, a_))
    colnames(xa_) <- colnames(df)
    muhat_ <- predict(mumod, newdata = xa_, type = "response")
    # gps_ <- dnorm(a_, x%*%beta[i-1,], sqrt(sigma2[i-1]))
    
    log.eps <- dpois(y, family$linkinv(family$linkfun(muhat_) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
      dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
      dpois(y, family$linkinv(family$linkfun(muhat[i - 1,]) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
      dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)
    
    test <- log(runif(n))
    a <- amat[i,] <- ifelse(((test <= log.eps) & !is.na(log.eps)), a_, a)
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
    xa <- cbind(x, predict(nsa, a))
    df <- data.frame(y = y, xa[,-1])

    xa.new.list <- lapply(a.vals, function(a.tmp, ...) {

      cbind(x, matrix(rep(c(predict(nsa, a.tmp)), n), byrow = T, nrow = n))

    })

    xa.new <- rbind(xa, do.call(rbind, xa.new.list))
    x.new <- as.matrix(xa.new[,1:ncol(x)])
    xa.new <- data.frame(xa.new[,-1])
    a.new <- c(a, rep(a.vals, each = n))
    colnames(x.new) <- colnames(x)
    colnames(xa.new) <- colnames(df)[-1]

    # exposure models
    pimod.vals <- c(x.new %*% beta[i,])
    pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2[i]))
    pihat <- pihat.vals[1:n]
    pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
    phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
    phat[which(phat < 0)] <- 1e-6

    # outcome models
    mumod <- glm(y ~ ., data = df, family = poisson, offset = offset)
    muhat.vals <- predict(mumod, newdata = xa.new, type = "response")
    muhat[i,] <- muhat.vals[1:n]
    muhat.mat <- matrix(muhat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
    mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y

    # integrate
    phat.mat <- matrix(rep(colMeans(pihat.mat), n), byrow = T, nrow = n)
    mhat.mat <- matrix(rep(colMeans(muhat.mat), n), byrow = T, nrow = n)
    intfn <- (muhat.mat - mhat.mat) * phat.mat
    int <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
                   (intfn[,-1] + intfn[,-length(a.vals)]) / 2, 1, sum)

    psi <- (y_ - muhat[i,])/(pihat/phat) + mhat

    dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = int, span = span, se.fit = TRUE)
    est.mat[i,] <- dr_out[1,]
    var.mat[i,] <- dr_out[2,]
    
    # xa <- cbind(1, predict(nsa, a))
    # df <- data.frame(y = y, xa)
    # gps <- dnorm(a, x%*%beta[i,], sqrt(sigma2[i]))
    # gamma_ <- gamma0 <- gamma[i-1,]
    
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
    # gamma[i,] <- gamma0
    
    # df <- data.frame(y = y, xa[,-1])
    # gmod <- rstanarm::stan_glm(y ~ ., data = df, family = poisson, QR = TRUE,
    #                            offset = offset, iter = 2*n.boot, chains = 1, refresh = 0)
    # gamma[i,,] <- do.call(cbind, split.along.dim(as.array(gmod), 3))
    
  }
  
  accept.a <- mean(apply(amat[(n.adapt + 1):nrow(amat),], 2, function(x) mean(diff(x) != 0) ))
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  gamma <- gamma[keep,,]
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  amat <- amat[keep,]
  est.mat <- est.mat[keep,]
  var.mat <- var.mat[keep,]
  
  # dr_out <- lapply(1:nrow(beta), function(i, ...){
  #   
  #   out <- sapply(a.vals, function(a.tmp, ...) {
  # 
  #     gps <- dnorm(a.tmp, x%*%beta[i,], sqrt(sigma2[i]))
  #     xa <- cbind(1, matrix(rep(predict(nsa, a.tmp), n), byrow = T, nrow = n), gps)
  #     val <- colMeans(family$linkinv(tcrossprod(xa, gamma[i,,])))
  #     return(list(est = mean(val), var = var(val)))
  # 
  #   })
  # 
  #   return(list(est = unlist(out[1,]), var = unlist(out[2,])))
  # 
  # })
  
  amat <- amat[,order(shield)]
  # est.mat <- do.call(rbind, lapply(dr_out, function(arg, ...) arg$est))
  # var.mat <- do.call(rbind, lapply(dr_out, function(arg, ...) arg$var))
  estimate <- colMeans(est.mat)
  variance <- colMeans(var.mat) + (1 + 1/nrow(amat))*apply(est.mat, 2, var)
  
  rslt <- list(estimate = estimate, variance = variance, accept.a = accept.a,
               mcmc = list(amat = amat, gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}

