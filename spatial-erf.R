
bart_spatial <- function(s, s.tilde, y, s.id, id, w = NULL, x = NULL, offset = NULL,
                         V = diag(1, length(y)), a.vals = seq(min(a), max(a), length.out = 100),
                         n.iter = 10000, n.adapt = 1000, thin = 10, 
                         shape = 1e-3, rate = 1e-3, scale = 1e6, 
                         bw = NULL, bw.seq = seq(0.1, 2, by = 0.1), folds = 5,
                         control = dbartsControl(updateState = FALSE, verbose = FALSE, n.burn = 0L, 
                                                 n.samples = 1L, n.trees = 200L, n.thin = 1L, n.chains = 1L)) {
  
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
  V <- V[shield, shield]
  
  if (is.null(w)) {
    ws <- cbind(rep(1, length(s.id)), s.tilde = s.tilde)
  } else {
    ws <- cbind(model.matrix(~ ., data.frame(w)), s.tilde = s.tilde)
  }
  
  s <- s[s.id %in% id]
  s.tilde <- s.tilde[s.id %in% id]
  ws <- w[s.id %in% id,,drop = FALSE]
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
  nu2 <- rep(NA, n.adapt + n.iter)
  rho <- rep(NA, n.adapt + n.iter)
  
  fit.a <- lm(a ~ 0 + x)
  fit.s <- lm(s.obs ~ 0 + ws.obs)
  
  beta[1,] <- coef(fit.a)
  alpha[1,] <- coef(fit.s)
  sigma2[1] <- sigma(fit.a)^2
  tau2[1] <- sigma(fit.s)^2
  omega2[1] <- var(s.hat - a.s)
  
  # initialize spatial stuff
  res.temp <- a - x %*% beta[1,]
  res.sd <- sd(res.temp, na.rm = TRUE)/5
  phi <- rnorm(n, mean = 0, sd = res.sd)
  nu2[1] <- var(phi) / 10
  rho[1] <- runif(1)
  accept.rho <- 0
  h.rho <- 0.02
  
  # CAR quantities
  V.quants <- check(V)
  V <- V.quants$W
  V3 <- V.quants$W.triplet
  n3 <- V.quants$n.triplet
  V3.sum <- V.quants$W.triplet.sum
  V.idx <- V.quants$W.begfin
  
  # determinant of V for updating rho
  Vstar <- diag(rowSums(V)) - V
  Vstar.eigen <- eigen(Vstar, symmetric = TRUE)
  Vstar.val <- Vstar.eigen$values
  detQ <- 0.5 * sum(log((rho[1] * Vstar.val + (1 - rho[1]))))
  
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
    
    print(i)
    
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
        dnorm(a_, phi + c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) +
        dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
        dnorm(ybar, samples$train, mean(samples$sigma)/sqrt(weights), log = TRUE) -
        dnorm(a, phi + c(x%*%beta[i - 1,]), sqrt(sigma2[i - 1]), log = TRUE) -
        dnorm(z.hat, a, sqrt(omega2[i - 1]/stab), log = TRUE)
      
      exam <- (log(runif(n)) <= log.eps) & !is.na(log.eps)
      temp <- ifelse(exam, a_, a)
      test <- sampler$setPredictor(x = temp, column = "a")
      
    }
    
    accept.a <- accept.a + exam
    a <- temp
    a.s <- rep(a, stab)
    
    # sample pred parameters
    alpha_var <- solve(t(ws.obs) %*% ws.obs + diag(tau2[i - 1]/scale, q, q))
    alpha[i,] <- rmvnorm(1, alpha_var %*% t(ws.obs) %*% s.obs, tau2[i - 1]*alpha_var)
    tau2[i] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                          sum(c(s.obs - c(ws.obs %*% alpha[i,]))^2)/2)
    
    # sample agg parameters
    omega2[i] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((s.hat - a.s)^2)/2)
    
    # sample gps parameters
    beta_var <- solve(t(x) %*% x + diag(sigma2[i - 1]/scale, p, p))
    beta[i,] <- rmvnorm(1, beta_var %*% t(x) %*% (a - phi), sigma2[i - 1]*beta_var)
    sigma2[i] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - phi - c(x %*% beta[i,]))^2)/2)
    
    # sample phi
    off.phi <- (a - as.numeric(x %*% beta[i,])) / sigma2[i]    
    phi <- CARBayes:::gaussiancarupdate(Wtriplet = V3, Wbegfin = V.idx, 
                                        Wtripletsum = V3.sum, nsites = n,
                                        phi = phi, tau2 = nu2[i - 1], rho = rho[i - 1], 
                                        nu2 = sigma2[i], offset = off.phi)
    phi <- phi - mean(phi)
    
    # sample nu2
    temp <- CARBayes:::quadform(V3, V3.sum, n3, n, phi, phi, rho[i - 1])
    nu2[i] <- 1/rgamma(1, shape + n/2, rate = rate + temp)
    
    # sample rho
    rho_ <- truncnorm::rtruncnorm(n = 1, a = 0, b = 1, mean = rho[i-1], sd = h.rho)  
    temp_ <- CARBayes:::quadform(V3, V3.sum, n3, n, phi, phi, rho_)
    detQ_ <- 0.5 * sum(log((rho_ * Vstar.val + (1 - rho_))))              
    logprob <- detQ - temp / nu2[i]
    logprob_ <- detQ_ - temp_ / nu2[i]
    hastings <- log(truncnorm::dtruncnorm(x = rho[i - 1], a = 0, b = 1, mean = rho_, sd = h.rho)) -
      log(truncnorm::dtruncnorm(x = rho_, a = 0, b = 1, mean = rho[i - 1], sd = h.rho)) 
    log.eps <- logprob_ - logprob + hastings
    
    if ((log(runif(1)) <= log.eps) & !is.na(log.eps)) {
      rho[i] <- rho_
      detQ <- detQ_
      accept.rho <- accept.rho + 1  
    } else
      rho[i] <- rho[i - 1]
    
    # sample outcome tree
    samples <- sampler$run()
    
    if (ceiling(i/100) == floor(i/100) & i < n.adapt) {
      
      h.rho <- ifelse(accept.rho > 40, h.rho + 0.1 * h.rho,
                      ifelse(accept.rho < 30, h.rho - 0.1 * h.rho, h.rho))
      accept.rho <- 0
      
      h.a <- ifelse(accept.a > 40, h.a + 0.1 * h.a,
                    ifelse(accept.a < 20, h.a - 0.1 * h.a, h.a))
      accept.a <- rep(0, n)
      
    }
    
    # save output
    if (i > n.adapt & (i - n.adapt)%%thin == 0) {
      
      j <- (i - n.adapt)/thin
      a.mat[j,] <- a 
      
      # outcome model
      muhat <- rowMeans(samples$train)
      
      muhat.mat <- sapply(a.vals, function(a.tmp, ...){
        
        xa.tmp <- data.frame(x[,-1], a = rep(a.tmp, n))
        colnames(xa.tmp) <- colnames(xa.train)[-1]
        c(sampler$predict(xa.tmp))
        
      })
      
      mhat.vals <- apply(muhat.mat, 2, weighted.mean, w = weights)
      mhat <- predict(smooth.spline(x = a.vals, y = mhat.vals), x = a)$y
      mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)
      
      # pseudo-outcome
      psi.mat[j,] <- psi <- (ybar - muhat + mhat)
      
      # integration matrix
      pimod.vals <- c(phi + x %*% beta[i,])
      
      # density estimation
      phat.vals <- sapply(a.vals, function(a.tmp, ...){
        mean(dnorm(a.tmp, pimod.vals, sqrt(sigma2[i])))
      })
      
      phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)  
      int.mat <- (muhat.mat - mhat.mat)*phat.mat
      
      # select bw if NULL
      if (j == 1 & is.null(bw))
        bw <- cv_bw(a = a, psi = psi, folds = folds, bw.seq = bw.seq)
      
      # asymptotics
      out <- sapply(a.vals, kern_est, psi = psi, a = a, weights = weights, 
                    bw = bw, se.fit = TRUE, int.mat = int.mat, a.vals = a.vals)
      
      est.mat[j,] <- out[1,]
      var.mat[j,] <- out[2,]
      
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
