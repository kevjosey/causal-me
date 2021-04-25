
gibbs_dr <- function(s, star, y, s.id, id, family = gaussian(),
                     offset = rep(0, length(id)), w = NULL, x = NULL,
                     shape = 1e-3, rate = 1e-3, scale = 1e6,
                     thin = 10, n.iter = 10000, n.adapt = 1000,
                     h.a = 0.5, h.gamma = 0.1, deg.num = 3, span = 0.75,
                     a.vals = seq(min(a), max(a), length.out = 20), mc.cores = 4) {
  
  # remove any s.id not present in id
  su.id <- unique(s.id)[order(unique(s.id))]
  su.id <- su.id[su.id %in% id]
  
  if(length(su.id) != length(id))
    stop("some observations in y.id are not represented by measurements of s.id.\n
         There are no exposure data for these entries.")
  
  if(!all(su.id == id[order(id)]))
    stop("some observations in y.id are not represented by measurements of s.id.\n
         There are no exposure data for these entries.")
  
  s <- s[s.id %in% id]
  star <- star[s.id %in% id]
  w <- w[s.id %in% id,]
  s.id <- s.id[s.id %in% id]
  
  m <- length(s.id)
  n <- length(id)
  
  if (is.null(x)) {
    x <- matrix(1, nrow = length(id), ncol = 1)
  } else {
    x <- model.matrix(~ ., data.frame(x))
  }
  
  if (is.null(w)) {
    ws <- cbind(rep(1, length(s.id)), ns(star, deg.num))
  } else {
    ws <- model.matrix(~ . + ns(star, deg.num), data.frame(w))
  }
  
  ws.tmp <- ws[!is.na(s),]
  s.tmp <- s[!is.na(s)]
  
  # initialize exposures
  a.tmp <- predict(lm(s.tmp ~ 0 + ., data = data.frame(ws.tmp)), newdata = data.frame(ws))
  a <- aggregate(a.tmp, by = list(s.id), mean)[,2]
  nsa <- ns(a, deg.num, Boundary.knots = c(min(a) - 2*sd(a)/sqrt(n), max(a) + 2*sd(a)/sqrt(n)))
  xa <- as.matrix(cbind(x, nsa))
  a.s <- rep(NA, length(s.id))
  
  stab <- table(s.id)
  sword <- order(order(s.id))
  shield <- order(order(id))
  a.s <- rep(a, stab)[sword]
  
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
  
  gamma[1,] <- coef(glm(y ~ 0 + xa, family = family, offset = offset))
  beta[1,] <- coef(lm(a ~ 0 + x))
  alpha[1,] <- coef(lm(s.tmp ~ 0 + ws.tmp))
  sigma2[1] <- sigma(lm(a ~ 0 + x))^2
  tau2[1] <- sigma(lm(s.tmp ~ 0 + ws.tmp))^2
  omega2[1] <- var(star - a.s)
  
  amat <- matrix(NA, nrow = n.iter + n.adapt, ncol = n)
  
  # gibbs sampler for predictors
  for(j in 2:(n.iter + n.adapt)) {
    
    # sample S
    
    sig <- sqrt((1/omega2[j - 1] + 1/tau2[j - 1])^(-1))
    hat <- (sig^2)*(a.s/omega2[j - 1] + c(ws %*% alpha[j - 1,])/tau2[j - 1])
    s.hat <- rnorm(m, hat, sig)
    s.hat[!is.na(s)] <- s[!is.na(s)]

    # sample A
    
    a_ <-  rnorm(n, a, h.a)
    xa_ <- as.matrix(cbind(x, predict(nsa, a_)))
    dt <- setDT(data.frame(s.id = s.id, s.hat = s.hat))[,.(A = .N, B = mean(s.hat)), by = 's.id'][order(s.id)]
    z.hat <- dt$B[shield]
    size <- dt$A[shield]
    
    log.eps <- dpois(y, exp(c(xa_%*%gamma[j - 1,]) + offset), log = TRUE) +
      dnorm(a_, c(x%*%beta[j - 1,]), sqrt(sigma2[j - 1]), log = TRUE) + 
      dnorm(a_, z.hat, sqrt(omega2[j - 1]/size), log = TRUE) -
      dpois(y, exp(c(xa%*%gamma[j - 1,]) + offset), log = TRUE) -
      dnorm(a, c(x%*%beta[j - 1,]), sqrt(sigma2[j - 1]), log = TRUE) -
      dnorm(a, z.hat, sqrt(omega2[j - 1]/size), log = TRUE)
    
    a <- amat[j,] <- ifelse((log(runif(1)) <= log.eps) & !is.na(log.eps), a_, a)
    
    xa <- as.matrix(cbind(x, predict(nsa, a)))
    a.s <- rep(a, stab)[sword]
    
    # Sample pred parameters
    
    alpha_var <- solve(t(ws.tmp) %*% ws.tmp + diag(tau2[j - 1]/scale, q, q))
    alpha[j,] <- rmvnorm(1, alpha_var %*% t(ws.tmp) %*% s.tmp, tau2[j - 1]*alpha_var)
    
    tau2[j] <- 1/rgamma(1, shape = shape + l/2, rate = rate +
                            sum(c(s.tmp - c(ws.tmp %*% alpha[j,]))^2)/2)
    
    # Sample agg parameters
    
    omega2[j] <- 1/rgamma(1, shape = shape + m/2, rate = rate + sum((s.hat - a.s)^2)/2)
    
    # Sample GPS parameters
    
    beta_var <- solve(t(x) %*% x + diag(sigma2[j - 1]/scale, p, p))
    beta[j,] <- rmvnorm(1, beta_var %*% t(x) %*% a, sigma2[j - 1]*beta_var)
    
    sigma2[j] <- 1/rgamma(1, shape = shape + n/2, rate = rate + sum(c(a - x %*% beta[j,])^2)/2)
   
    gamma_ <- gamma[j,] <- gamma[j-1,]
    
    for (k in 1:ncol(xa)){
      
      gamma_[k] <- c(rnorm(1, mean = gamma[j,k], h.gamma))
      
      log.eps <- sum(dpois(y, exp(c(xa %*% gamma_) + offset), log = TRUE)) -
        sum(dpois(y, exp(c(xa %*% gamma[j,]) + offset), log = TRUE)) +
        dnorm(gamma_[k], 0, scale, log = TRUE) - dnorm(gamma[j,k], 0 , scale, log = TRUE)
        
      if ((log(runif(1)) <= log.eps) & !is.na(log.eps))
        gamma[j,k] <- gamma_[k]
      else
        gamma_[k] <- gamma[j,k]
      
    }
    
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
  
  y.new <- exp(log(y) - offset)

  out <- mclapply(1:nrow(amat), function(k, ...){

    a <- amat[k,]
    xa <- as.matrix(cbind(x, predict(nsa, a)))
    
    if (ncol(x) > 1) {
      pimod.vals <- c(x %*% beta[k,])
    } else {
      pimod.vals <- c(x %*% beta[k])
    }

    # exposure models
    pihat <- dnorm(a, pimod.vals, sqrt(sigma2[k]))
    phat.vals <- sapply(a.vals, function(a.tmp, ...)
      mean(dnorm(a.tmp, pimod.vals, sqrt(sigma2[k]))))
    phat <- predict(smooth.spline(a.vals, phat.vals), x = a)$y
    phat[which(phat < 0)] <- 1e-6
    phat.mat <- matrix(rep(phat.vals, n), byrow = T, nrow = n)

    muhat <- exp(c(xa %*% gamma[k,]))

    # predict marginal outcomes given a.vals (or a.agg)
    muhat.mat <- sapply(a.vals, function(a.tmp, ...) {

      xa.tmp <- cbind(x = x, matrix(rep(c(predict(nsa, a.tmp)), n), 
                                    byrow = TRUE, nrow = n))
      return(exp(c(xa.tmp %*% gamma[k,])))

    })

    # aggregate muhat.vals and integrate for influence function
    mhat.vals <- colMeans(muhat.mat)
    mhat <- predict(smooth.spline(a.vals, mhat.vals), x = a)$y
    mhat.mat <- matrix(rep(mhat.vals, n), byrow = T, nrow = n)

    # integrate
    intfn <- (muhat.mat - mhat.mat) * phat.mat
    int <- apply(matrix(rep((a.vals[-1]-a.vals[-length(a.vals)]), n), byrow = T, nrow = n) *
                   (intfn[,-1] + intfn[,-length(a.vals)]) / 2, 1, sum)

    psi <- c((y.new - muhat)/(pihat/phat) + mhat)
    
    dr_out <- sapply(a.vals, dr_est, psi = psi, a = a, int = int,
                     family = family, span = span, se.fit = FALSE)

    return(dr_out)
           
  }, mc.cores = mc.cores)

  est.mat <- do.call(rbind, out)
  estimate <- colMeans(est.mat)
  variance <- apply(est.mat, 2, var)
  
  amat <- amat[,order(shield)]

  rslt <- list(estimate = estimate, variance = variance,
               mcmc = list(gamma = gamma, beta = beta, alpha = alpha,
                           sigma2 = sigma2, tau2 = tau2, omega2 = omega2,
                           amat = amat),
               accept.a = accept.a, accept.gamma = accept.gamma)
  
  return(rslt)
  
}
