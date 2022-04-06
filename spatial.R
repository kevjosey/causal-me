
bart_spatial <- function(s, t, y, s.id, id, w = NULL, x = NULL, V = diag(1, length(y)),
                         family = gaussian(), offset = NULL, weights = NULL,
                         a.vals = seq(min(a), max(a), length.out = 100),
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
  
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
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
  h.a <- 0.5*sd(a)
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
  
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(s.hat - a.s)
  
  # initialize spatial stuff
  res.temp <- y - x %*% beta[1,]
  res.sd <- sd(res.temp, na.rm = TRUE)/5
  phi <- rnorm(n, mean = 0, sd = res.sd)
  nu2[i] <- var(phi) / 10
  rho[i] <- runif(1)
  accept.rho <- c(0,0)
  h.rho <- 0.02
  
  # CAR quantities
  V.quants <- CARBayes:::common.Wcheckformat(V)
  V <- V.quants$W
  V3 <- V.quants$W.triplet
  n3 <- V.quants$n.triplet
  V3.sum <- V.quants$W.triplet.sum
  V.idx <- V.quants$W.begfin
  
  # determinant of V for updating rho
  Vstar <- diag(apply(V, 1, sum)) - V
  Vstar.eigen <- eigen(Vstar)
  Vstar.val <- Vstar.eigen$values
  detQ <- 0.5 * sum(log((rho[1] * Vstar.val + (1 - rho[1]))))
  
  # initialize bart
  ybar <- family$linkinv(family$linkfun(y) - offset)
  xa.train <- data.frame(ybar = ybar, x[,-1], a = a)
  sampler <- dbarts::dbarts(ybar ~ ., data = xa.train, control = control, weights = weights)
  
  # run first iteration of tree
  samples <- sampler$run()
  
  # ERC stuff
  a.mat <- psi <- matrix(NA, nrow = floor(n.iter/thin), ncol = n)
  mhat.out <- est.mat <- var.mat <- matrix(NA, nrow = floor(n.iter/thin), ncol = length(a.vals))
  
  # gibbs sampler for predictors
  for(i in 2:(n.iter + n.adapt)) {
    
    # print(i)
    
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
        dnorm(a_, z.hat, sqrt(omega2[i - 1]/stab), log = TRUE) -
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
    alpha[i,] <- rmvnorm(1, alpha_var %*% t(ws.tmp) %*% s.tmp, tau2[i - 1]*alpha_var)
    tau2[i] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                          sum(c(s.tmp - c(ws.tmp %*% alpha[i,]))^2)/2)
    
    # sample agg parameters
    omega2[i] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((s.hat - a.s)^2)/2)
    
    # sample gps parameters
    beta_var <- solve(t(x) %*% x + diag(sigma2[i - 1]/scale, p, p))
    beta[i,] <- rmvnorm(1, beta_var %*% t(x) %*% (a - phi), sigma2[i - 1]*beta_var)
    sigma2[i] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - phi - c(x %*% beta[i,]))^2)/2)
    
    # sample phi
    off.phi <- (ybar - as.numeric(x %*% beta[i,])) / sigma2[i]    
    phi <- CARBayes:::gaussiancarupdate(Wtriplet = V3, Wbegfin = V.idx, 
                                        Wtripletsum = V3.sum, nsites = n,
                                        phi = phi, tau2 = nu2[i], rho = rho[i], 
                                        nu2 = sigma2[i], offset = off.phi)
    phi <- phi - mean(phi)
    
    # sample nu2
    temp <- CARBayes:::quadform(V3, V3.sum, n3, n, phi, phi, rho[i])
    nu2[i] <- 1/rgamma(1, shape + n/2, rate = rate + temp)
    
    # sample rho
    rho_ <- truncnorm::rtruncnorm(n = 1, a = 0, b = 1, mean = rho[i], sd = h.rho)  
    temp_ <- CARBayes:::quadform(V3, V3.sum, n3, n, phi, phi, rho_)
    detQ_ <- 0.5 * sum(log((rho_ * Vstar.val + (1 - rho_))))              
    logprob <- detQ - temp / nu2[i]
    logprob_ <- detQ_ - temp_ / nu2[i]
    hastings <- log(truncnorm::dtruncnorm(x = rho[i - 1], a = 0, b = 1, mean = rho_, sd = h.rho)) -
      log(truncnorm::dtruncnorm(x = rho_, a = 0, b = 1, mean = rho[i - 1], sd = h.rho)) 
    log.eps <- logprob_ - logprob + hastings
    accept.rho[2] <- accept.rho[2] + 1  
    
    if ((log(runif(1)) <= log.eps) & !is.na(log.eps)) {
      rho[i] <- rho_
      detQ <- detQ_
      accept.rho[1] <- accept.rho[1] + 1  
    } else
      rho[i] <- rho[i - 1]

    # sample outcome tree
    samples <- sampler$run()
    
    if (ceiling(i/100) == floor(i/100) & j < n.adapt) {
      
      h.rho <- CARBayes:::common.accceptrates2(accept.rho, h.rho, 30, 40, 0.5)
      accept.rho <- c(0,0)
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
      
      muhat.mat <- sapply(a.vals, function(a.tmp, ...){
        xa.tmp <- data.frame(x[,-1], a = rep(a.tmp, n))
        colnames(xa.tmp) <- colnames(xa.train)[-1]
        c(sampler$predict(xa.tmp))
      })
      
      mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
      
      # pseudo-outcome
      psi[j,] <- (ybar - muhat + mhat)
      mhat.out[j,] <- colMeans(muhat.mat)
      
      # exposure model for integtion
      a.std <- c(c(a, a.vals) - mean(a)) / sd(a)
      dens <- density(a.std[1:n])
      phat.vals <- approx(x = dens$x, y = dens$y, xout = a.std[-(1:n)])$y / sd(a)
      
      # integration matrix
      phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)
      mhat.mat <- matrix(rep(colMeans(muhat.mat), n), byrow = T, nrow = n)
      int.mat <- (muhat.mat - mhat.mat) * phat.mat
      
      if (j == 1 & is.null(span))
        span <- cv_span(a = a, psi = psi[j,], family = family, folds = folds, span.seq = span.seq)
      
      out <- sapply(a.vals, kern_est, psi = psi[j,], a = a, span = span, 
                    family = gaussian(), se.fit = TRUE, int.mat = int.mat, 
                    a.vals = a.vals)
      
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
  smooth_estimate <- colMeans(est.mat)
  smooth_variance <- colMeans(var.mat) + (1 + 1/nrow(a.mat))*apply(est.mat, 2, var)
  tree_estimate <- colMeans(mhat.out)
  tree_variance <- apply(mhat.out, 2, var)
  hpdi <- apply(mhat.out, 2, hpd)
  rownames(hpdi) <- c("lower", "upper")
  
  rslt <- list(smooth_estimate = smooth_estimate, smooth_variance = smooth_variance,
               tree_estimate = tree_estimate, tree_variance = tree_variance, hpdi = hpdi,
               accept.a = accept.a, mcmc = list(a.mat = a.mat, beta = beta, alpha = alpha, 
                                                sigma2 = sigma2, tau2 = tau2, omega2 = omega2))
  
  return(rslt)
  
}
