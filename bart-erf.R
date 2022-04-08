
bart_erf <- function(s, t, y, s.id, id, w = NULL, x = NULL, offset = NULL,
                     a.vals = seq(min(s), max(s), length.out = 100),
                     n.iter = 10000, n.adapt = 1000, thin = 10, 
                     shape = 1e-3, rate = 1e-3, scale = 1e6,
                     span = NULL, span.seq = seq(0.1, 2, by = 0.1), folds = 5,
                     control = dbartsControl(updateState = FALSE, verbose = FALSE, n.burn = 0L, 
                                             n.samples = 1L, n.thin = thin, n.chains = 1L)) {
  
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
  
  s <- s[s.id %in% id]
  t <- as.matrix(t)[s.id %in% id]
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
  offset <- offset[shield]
  weights <- weights[shield]
  
  if (is.null(w)) {
    ws <- cbind(rep(1, length(s.id)), t = t)
  } else {
    ws <- cbind(model.matrix(~ ., data.frame(w)), t = t)
  }
  
  sword <- order(s.id)
  stab <- table(s.id)
  s <- s[sword]
  t <- t[sword]
  ws <- ws[sword,]
  s.id <- s.id[sword]
  
  # when s is observed
  ws.tmp <- ws[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  
  # initialize exposures
  s.hat <- predict(lm(s.tmp ~ 0 + ., data = data.frame(ws.tmp)), newdata = data.frame(ws))
  a <- aggregate(s.hat, by = list(s.id), mean)[,2]
  a.s <- rep(a, stab)
  h.a <- rep(0.5*sd(a), n)
  accept.a <- rep(0, n)

  # data dimensions
  p <- ncol(x)
  q <- ncol(ws)
  
  # initialize parameters
  beta <- matrix(NA, nrow = n.adapt + n.iter, ncol = p)
  alpha <- matrix(NA, nrow = n.adapt + n.iter, ncol = q)
  sigma2 <- rep(NA, n.adapt + n.iter)
  tau2 <- rep(NA, n.adapt + n.iter)
  omega2 <- rep(NA, n.adapt + n.iter)
  
  fit.a <- lm(a ~ 0 + x)
  fit.s <- lm(s.tmp ~ 0 + ws.tmp)
  
  beta[1,] <- coef(fit.a)
  alpha[1,] <- coef(fit.s)
  sigma2[1] <- sigma(fit.a)^2
  tau2[1] <- sigma(fit.s)^2
  omega2[1] <- var(s.hat - a.s)
  
  # initialize bart
  ybar <- exp(log(y) - offset)
  xa.train <- data.frame(ybar = ybar, x[,-1], a = a)
  sampler <- dbarts::dbarts(ybar ~ ., data = xa.train, control = control, weights = weights)
  
  # run first iteration of tree
  samples <- sampler$run()
  
  # initialize output
  a.mat <- psi.mat <- matrix(NA, nrow = floor(n.iter/thin), ncol = n)
  est.mat <- var.mat <- matrix(NA, nrow = floor(n.iter/thin), ncol = length(a.vals))
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    # sample S
    sig.s <- sqrt((1/omega2[i - 1] + 1/tau2[i - 1])^(-1))
    mu.s <- (sig.s^2)*(a.s/omega2[i - 1] + c(ws %*% alpha[i - 1,])/tau2[i - 1])
    s.hat <- rnorm(m, mu.s, sig.s)
    s.hat[!is.na(s)] <- s[!is.na(s)]
    
    # sample A and update outcome tree
    test <- FALSE
    z.hat <- aggregate(s.hat, by = list(s.id), mean)[,2]
    
    while (test == FALSE) {
      
      a_ <- rnorm(n, a, h.a)
      xa.test <- data.frame(x = x[,-1], a = a_)
      colnames(xa.test) <- colnames(xa.train)[-1]
      xa.pred <- sampler$predict(xa.test)
      
      log.eps <- dnorm(ybar, xa.pred, mean(samples$sigma)/sqrt(weights), log = TRUE) +
        dnorm(a_, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
        dnorm(z.hat, a_, sqrt(omega2[i - 1]/stab), log = TRUE) -
        dnorm(ybar, samples$train, mean(samples$sigma)/sqrt(weights), log = TRUE) -
        dnorm(a, c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
        dnorm(z.hat, a, sqrt(omega2[i - 1]/stab), log = TRUE)
      
      exam <- (log(runif(n)) <= log.eps) & !is.na(log.eps)
      temp <- ifelse(exam, a_, a)
      test <- sampler$setPredictor(x = temp, column = "a")
      
    }
    
    accept.a <- accept.a + exam
    a <- temp
    a.s <- rep(a, stab)
    
    # sample pred parameters
    alpha_var <- solve(t(ws.tmp) %*% ws.tmp + diag(tau2[i - 1]/scale, q, q))
    alpha[i,] <- rmvnorm(1, alpha_var %*% c(t(ws.tmp) %*% s.tmp), tau2[i - 1]*alpha_var)
    tau2[i] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                          sum(c(s.tmp - c(ws.tmp %*% alpha[i,]))^2)/2)
    
    # sample agg parameters
    omega2[i] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((s.hat - a.s)^2)/2)
    
    # sample gps parameters
    beta_var <- solve(t(x) %*% x + diag(sigma2[i - 1]/scale, p, p))
    beta[i,] <- rmvnorm(1, beta_var %*% c(t(x) %*% a), sigma2[i - 1]*beta_var)
    sigma2[i] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - c(x %*% beta[i,]))^2)/2)
    
    # sample outcome tree
    samples <- sampler$run()
    
    # update random walk parameters
    if(ceiling(i/100) == floor(i/100) & i < n.adapt) {
      
      h.a <- ifelse(accept.a > 20, h.a + 0.1 * h.a,
                    ifelse(accept.a < 10, h.a - 0.1 * h.a, h.a))
      accept.a <- rep(0, n)
      
    }
    
    # save output
    if (i > n.adapt & (i - n.adapt)%%thin == 0) {
      
      j <- (i - n.adapt)/thin
      a.mat[j,] <- a 
      
      # outcome model
      muhat <- rowMeans(samples$train)
      
      muhat.mat <- sapply(a, function(a.tmp, ...){
        
        xa.tmp <- data.frame(x[,-1], a = rep(a.tmp, n))
        colnames(xa.tmp) <- colnames(xa.train)[-1]
        c(sampler$predict(xa.tmp))
        
      })

      mhat <- apply(muhat.mat, 2, weighted.mean, w = weights)

      # pseudo-outcome
      psi.mat[j,] <- psi <- (ybar - muhat + mhat)
      
      # integration matrix
      int.mat <- muhat.mat - matrix(rep(mhat, n), byrow = T, nrow = n)

      # select span if NULL
      if (j == 1 & is.null(span))
        span <- cv_span(a = a, psi = psi, folds = folds, span.seq = span.seq)

      # asymptotics
      out <- sapply(a.vals, loess_est, psi = psi, a = a, weights = weights, 
                    span = span, se.fit = TRUE, int.mat = int.mat)

      est.mat[j,] <- out[1,]
      var.mat[j,] <- out[2,]
      
      # bootstrap
      # boot.idx <- cbind(1:n, replicate(200, sample(x = n, size = n, replace = TRUE)))
      # 
      # out <- apply(boot.idx, 2, function(idx, ...) {
      # 
      #   sapply(a.vals, loess_est, psi = psi[idx], a = a[idx], 
      #          weights = weights[idx], span = span, se.fit = FALSE)
      # 
      # })
      # 
      # est.mat[j,] <- out[,1]
      # var.mat[j,] <- apply(out[,2:ncol(out)], 1, var)
      
    }
    
  }
  
  accept.a <- apply(a.mat, 2, function(x) mean(diff(x) != 0) )
  
  # thinning
  keep <- seq(n.adapt + 1, n.iter + n.adapt, by = thin)
  beta <- beta[keep,]
  alpha <- alpha[keep,]
  sigma2 <- sigma2[keep]
  tau2 <- tau2[keep]
  omega2 <- omega2[keep]
  
  a.mat <- a.mat[,order(shield)]
  psi.mat <- psi.mat[,order(shield)]
  estimate <- colMeans(est.mat)
  variance <- colMeans(var.mat) + (1 + 1/nrow(a.mat))*apply(est.mat, 2, var)
  
  rslt <- list(estimate = estimate, variance = variance, accept.a = accept.a,
               mcmc = list(a.mat = a.mat, psi.mat = psi.mat, beta = beta, alpha = alpha, 
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}