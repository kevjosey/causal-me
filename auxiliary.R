
ipw <- function(a, x, beta, sigma2, a.vals) {
  
  n <- length(a)
  x.new <- x[rep(1:n, length(a.vals) + 1), ]
  a.new <- c(a, rep(a.vals, each = n))
  pimod.vals <- c(x.new %*% beta)
  pihat.vals <- dnorm(a.new, pimod.vals, sqrt(sigma2))
  pihat <- pihat.vals[1:n]
  pihat.mat <- matrix(pihat.vals[-(1:n)], nrow = n, ncol = length(a.vals))
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat)), x = a)$y
  phat[which(phat < 0)] <- 1e-6
  out <- phat/pihat
  return(out)
  
}

split.along.dim <- function(a, n) {
  
  setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[,n]),
                  array, dim = dim(a)[-n], dimnames(a)[-n]),
           dimnames(a)[[n]])
  
}

# highest posterior density
hpd <- function(x, alpha = 0.05){
  
  n <- length(x)
  m <- round(n * alpha)
  x <- sort(x)
  y <- x[(n - m + 1):n] - x[1:m]
  z <- min(y)
  k <- which(y == z)[1]
  c(x[k], x[n - m + k])
  
}

mclapply.hack <- function(..., mc.cores = 1) {
  
  ## Create a cluster
  ## ... How many workers do you need?
  ## ... N.B. list(...)[[1]] returns the first
  ##          argument passed to the function. In
  ##          this case it is the list to iterate over
  size.of.list <- length(list(...)[[1]])
  cl <- makeCluster( mc.cores )
  
  ## Find out the names of the loaded packages
  loaded.package.names <- c(
    ## Base packages
    sessionInfo()$basePkgs,
    ## Additional packages
    names( sessionInfo()$otherPkgs ))
  
  ## N.B. tryCatch() allows us to properly shut down the
  ##      cluster if an error in our code halts execution
  ##      of the function. For details see: help(tryCatch)
  tryCatch( {
    
    ## Copy over all of the objects within scope to
    ## all clusters.
    ##
    ## The approach is as follows: Beginning with the
    ## current environment, copy over all objects within
    ## the environment to all clusters, and then repeat
    ## the process with the parent environment.
    ##
    this.env <- environment()
    while( identical( this.env, globalenv() ) == FALSE ) {
      clusterExport(cl,
                    ls(all.names=TRUE, env=this.env),
                    envir=this.env)
      this.env <- parent.env(environment())
    }
    ## repeat for the global environment
    clusterExport(cl, ls(all.names=TRUE, env  = globalenv()), envir = parent.frame())
    
    ## Load the libraries on all the clusters
    ## N.B. length(cl) returns the number of clusters
    parLapply( cl, 1:length(cl), function(xx){
      lapply(loaded.package.names, function(yy) {
        ## N.B. the character.only option of
        ##      require() allows you to give the
        ##      name of a package as a string.
        require(yy , character.only=TRUE)})
    })
    
    ## Run the lapply in parallel
    return( parLapply( cl, ...) )
  }, finally = {
    ## Stop the cluster
    stopCluster(cl)
  })
}

hush <- function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

message_parallel <- function(...){
  system(sprintf('echo "\n%s\n"', paste(..., collapse="")))
}
