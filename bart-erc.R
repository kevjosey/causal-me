
bart_erc <- function(s, star, y, s.id, id, family = gaussian(),
                    offset = rep(0, length(id)), w = NULL, x = NULL,
                    a.vals = seq(min(a), max(a), length.out = 100),
                    shape = 1e-3, rate = 1e-3, scale = 1e6,
                    thin = 10, n.iter = 10000, n.adapt = 1000,
                    h.a = 0.5, span = 0.75, 
                    control = dbartsControl(updateState = FALSE, verbose = FALSE, n.burn = 0L, 
                                            n.samples = 1L, n.thin = thin, n.chains = 1L)) {
  
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
  a.mat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  
  # the good stuff
  a.mat <- int <- psi <- matrix(NA, nrow = n.iter, ncol = n)
  
  # initialize bart
  y_ <- family$linkinv(family$linkfun(y) - offset)
  xa.train <- data.frame(y_ = y_, x = x[,-1], a = a)
  xa.test <- data.frame(x = x[rep(1:n, length(a.vals) + 1),-1], a = c(a, rep(a.vals, each = n)))
  sampler <- dbarts::dbarts(y_ ~ ., xa.train, xa.test, control = control)
  
  # run first iteration of tree
  a.mat[1,] <- a
  a_ <- c(rnorm(n, a, h.a), rep(a.vals, each = n))
  sampler$setTestPredictor(x = a_, column = "a")
  samples <- sampler$run()
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    # print(i)
    
    # sample S
    
    sig.s <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    mu.s <- (sig.s^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, mu.s, sig.s)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A
    
    test <- FALSE
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    
    while (test == FALSE) {
      
      a_ <- c(rnorm(n, a, h.a), rep(a.vals, each = n))
      sampler$setTestPredictor(x = a_, column = "a")
      samples <- sampler$run()
      
      log.eps <- dnorm(y_, rowMeans(samples$test[1:n,,drop = FALSE]), mean(samples$sigma), log = TRUE) +
        dnorm(a_[1:n], c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
        dnorm(a_[1:n], z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
        dnorm(y_, rowMeans(samples$train), mean(samples$sigma), log = TRUE) -
        dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
        dnorm(a, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE)
      
      temp <- ifelse(((log(runif(n)) <= log.eps) & !is.na(log.eps)), a_[1:n], a)
      test <- sampler$setPredictor(x = temp, column = "a")
      
    }
    
    a <- temp
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
    
    # Sample outcome trees
    

    
    if (i > n.adapt) {
    
      offset <- i - n.adapt
      a.mat[offset,] <- a 
      
      # outcome model
      muhat <- rowMeans(samples$train)
      muhat.mat <- matrix(rowMeans(samples$test[-(1:n),,drop = FALSE]), nrow = n, ncol = length(a.vals))
      mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
      
      # exposure model
      pimod.vals <- c(x[rep(1:n, length(a.vals) + 1),] %*% beta[i,])
      pihat.vals <- dnorm(c(a, rep(a.vals, each = n)), pimod.vals, sqrt(sigma2[i]))
      pihat <- pihat.vals[1:n]
      pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
      phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
      phat[which(phat < 0)] <- 1e-6
      
      phat.mat <- matrix(rep(colMeans(pihat.mat), n), byrow = T, nrow = n)
      mhat.mat <- matrix(rep(colMeans(muhat.mat), n), byrow = T, nrow = n)
      intfn <- (muhat.mat - mhat.mat) * phat.mat
      int[offset,] <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
                     (intfn[,-1] + intfn[,-length(a.vals)]) / 2, 1, sum)
      
      psi[offset,] <- (y_ - muhat)/(pihat/phat) + mhat
      
    }
  
  }
  
  accept.a <- apply(a.mat[(n.adapt + 1):nrow(a.mat),], 2, function(x) mean(diff(x) != 0) )
  
  # thinning
  keep1 <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep1,]
  alpha <- alpha[keep1,]
  sigma2 <- sigma2[keep1]
  tau2 <- tau2[keep1]
  omega2 <- omega2[keep1]
  
  # thinning ERC input
  keep2 <- seq(1, n.iter, by = thin)
  a.mat <- a.mat[keep2,]
  psi <- psi[keep2,]
  int <- int[keep2,]
  
  # analyze output
  
  out <- lapply(1:nrow(a.mat), function(i, ...){
    
    a <- a.mat[i,]
    psi <- psi[i,]
    int <- int[i,]
    
    dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = int, span = span, se.fit = TRUE)
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
  
  rslt <- list(estimate = estimate, variance = variance, hpdi = hpdi, accept.a = accept.a,
               mcmc = list(a.mat = a.mat, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}
