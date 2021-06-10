
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
  
  # initialize parameters
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  a.mat <- z.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  est.mat <- var.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = length(a.vals))
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  a.mat[1,] <- z.mat[1,] <- a
  
  # outcome stuff
  y_ <- family$linkinv(family$linkfun(y) - offset)
  pihat <- dnorm(a, c(x %*% beta[1,]), sqrt(sigma2[1]), log = TRUE)
  xa <- cbind(a, pihat)
  mumod <- hush(kernlab::gausspr(y = y_, x = xa, kernel = "rbfdot", variance.model = TRUE))
  muhat <- predict(mumod, xa, type = "response")
  varhat <- predict(mumod, xa, type = "variance")
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    # print(i)
    
    # sample S
    
    sig.s <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    mu.s <- (sig.s^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, mu.s, sig.s)
    
    # sample A
    
    z.hat <- z.mat[i,] <- aggregate(s.hat, by = list(s.id), sum)[,2]
    sig.a <- sqrt((1/sigma2[i - 1] + stab/omega2[i - 1])^(-1))
    mu.a <- (sig.a^2)*(c(x %*% beta[i - 1,])/sigma2[i - 1] + z.hat/omega2[i - 1])
    a <- a.mat[i,] <- rnorm(n, mu.a, sig.a)
    a.s <- rep(a, stab)
    
    # z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    # a_ <- rnorm(n, a, h.a)
    # pihat_ <- dnorm(a_, c(x %*% beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE)
    # xa_ <- cbind(a_, pihat_)
    # colnames(xa_) <- colnames(xa)
    # muhat_ <- predict(mumod, xa_, type = "response")
    # varhat_ <- predict(mumod, xa_, type = "variance")
    # 
    # log.eps <- dnorm(y_, muhat_, sqrt(varhat_), log = TRUE) +
    #   dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
    #   dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
    #   dnorm(y_, muhat, sqrt(varhat), log = TRUE) -
    #   dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
    #   dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)
    # 
    # test <- log(runif(n))
    # a <- a.mat[i,] <- ifelse(((test <= log.eps) & !is.na(log.eps)), a_, a)
    # a.s <- rep(a, stab)
    
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
    
    sigma2[i] <- 1/rgamma(1, shape = shape + n/2, rate = rate + 
                            sum(c(a - c(x %*% beta[i,]))^2)/2)
    
    # Sample outcome model while cutting feedback
    
    # Metropolis-Hastings
    # pihat <- dnorm(a, c(x %*% beta[i,]), sqrt(sigma2[i]), log = TRUE)
    # xa <- cbind(a, pihat)
    # mumod <- hush(kernlab::gausspr(y = y_, x = xa, kernel = "rbfdot", variance.model = TRUE))
    # muhat <- predict(mumod, xa, type = "response")
    # varhat <- predict(mumod, xa, type = "variance")
    
  }
  
  accept.a <- mean(apply(a.mat[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) ))
  # accept.gamma <- mean(apply(gamma[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) ))
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  # gamma <- gamma[keep,,]
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
    sig.a <- sqrt((1/sigma2[i] + rep(stab, length(a.vals))/omega2[i])^(-1))
    mu.a <- (sig.a^2)*(c(x.new %*% beta[i,])/sigma2[i] + rep(z.mat[i,], length(a.vals))/omega2[i])
    # pimod.vals <- c(x.new %*% beta[i,])
    pihat.vals <- dnorm(a.new, mu.a, sig.a)
    pihat.mat <- matrix(pihat.vals, nrow = n, ncol = length(a.vals))
    phat <- colMeans(pihat.mat)

    # uncertainty design
    dr_out <- sapply(1:length(a.vals), function(j, ...){

      a.std <- a - a.vals[j]
      k <- floor(min(span, 1)*length(a))
      idx <- order(abs(a.std))[1:k]
      a.std <- a.std[idx]
      y <- y[idx]
      k.std <- pihat.mat[idx,j]/phat[j]
      gh <- cbind(1, a.std)
      b <- optim(par = c(0,0), fn = opt_fun2, k.std = k.std, psi = y, gh = gh, family = family, offset = offset[idx])
      mu <- family$linkinv(c(b$par[1]))

      gh.inv <- solve(t(gh) %*% diag(k.std) %*% gh)
      v.inf <- (y - family$linkinv(c(gh%*%b$par + offset[idx])))^2
      sig <- gh.inv %*% t(gh) %*% diag(k.std) %*% diag(v.inf) %*% diag(k.std) %*% gh %*% gh.inv
      return(c(mu = mu, sig = sig[1,1]))

    })

    estimate <- dr_out[1,]
    variance <- dr_out[2,]

    return(list(estimate = estimate, variance = variance))

  })
  
  # out <- lapply(1:nrow(a.mat), function(i, ...){
  # 
  #   a <- a.mat[i,]
  #   n <- length(a)
  # 
  #   # sig.a <- sqrt((1/sigma2[i] + rep(stab, length(a.vals))/omega2[i])^(-1))
  #   # mu.a <- (sig.a^2)*(c(x.new %*% beta[i,])/sigma2[i] + rep(z.mat[i,], length(a.vals))/omega2[i])
  # 
  #   pimod <- c(x %*% beta[i,])
  #   pihat <- dnorm(a, pimod, sqrt(sigma2[i]), log = TRUE)
  #   xa <- cbind(a, pihat)
  #   mumod <- hush(kernlab::gausspr(y = y_, x = xa, kernel = "rbfdot", variance.model = TRUE))
  #   
  #   dr_out <- sapply(a.vals, function(a.tmp, ...){
  #     
  #     pihat.vals <- dnorm(a.tmp, pimod, sqrt(sigma2[i]), log = TRUE)
  #     xa_ <- cbind(a.tmp, pihat.vals)
  #     muhat.vals <- predict(mumod, xa_, type = "response")
  #     var.vals <- predict(mumod, xa_, type = "variance")
  #     
  #     return(c(muhat.vals = mean(muhat.vals), var.vals = mean(var.vals)))
  #     
  #   })
  # 
  #   estimate <- dr_out[1,]
  #   variance <- dr_out[2,]
  # 
  #   return(list(estimate = estimate, variance = variance))
  # 
  # })
  
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
