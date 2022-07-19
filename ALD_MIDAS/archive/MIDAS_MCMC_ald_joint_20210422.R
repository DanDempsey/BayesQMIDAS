###### MCMC Routines for IRTS MIDAS Inference 
###### Daniel Dempsey

MIDAS_MCMC_ald <- function(formula, start, prior, quantile, MCMC_length) {
  
  ### Prepare parameters
  Zenv <- environment( formula )
  yv <- Zenv$response_vector
  n_y <- length(yv)
  rtrunc <- ifelse(yv, TRUE, FALSE)
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  lambda <- 2 + psi^2 / omega
  V0i <- chol2inv( chol(prior$V0) )
  V0ib0 <- V0i %*% prior$beta0
  
  X_mat <- model.matrix( formula, Zenv )
  nu <- rep(1, n_y)
  pn <- psi * nu
  V <- omega * nu
  
  beta <- matrix(0, ncol = ncol(X_mat), nrow = MCMC_length + 1L)
  colnames(beta) <- paste0( "beta_", 0L:(ncol(beta)-1L) )
  beta[1L, ] <- start

  ### Main Loop
  for ( i in 2L:nrow(beta) ) {
    
    ### Update z
    Xb <- X_mat %*% beta[i-1L, ]
    #z <- rtruncnorm( n = n_y, a = lower, b = upper, mean = Xb + pn, sd = sqrt(V) )
    z <- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
    
    ### Update nu
    mu <- sqrt( lambda * omega / (z - Xb)^2  )
    nu <- 1 / rinvgauss( n = n_y, mean = mu, shape = lambda )
    
    ### Update beta
    V <- omega * nu
    pn <- psi * nu
    XtVi <- t( X_mat / V )
    V_post <- chol2inv( chol(V0i + XtVi%*%X_mat) )
    B_post <- V_post%*%( V0ib0 + XtVi%*%(z - pn) )
    beta[i, ] <- mvrnorm( mu = B_post, Sigma = V_post )
    
    ### Print progress
    if ( !i%%500L )
      print( paste0("Current iteration: ", i) )

  }
  
  beta
  
}

MIDAS_MCMC_ald_vs <- function(formula, start, prior, quantile, MCMC_length) {
  
  ### Preparation
  Zenv <- environment( formula )
  yv <- Zenv$response_vector
  n_y <- length(yv)
  rtrunc <- ifelse(yv, TRUE, FALSE)
  
  X_mat_full <- model.matrix( formula, Zenv )
  V0i_full <- chol2inv( chol(prior$V0) )
  V0ib0_full <- V0i_full %*% prior$beta0
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  lambda <- 2 + psi^2 / omega
  nu <- rep(1, n_y)
  z <- rep(0, n_y)
  pn <- psi * nu
  V <- omega * nu
  vs_prob <- prior$vars
  n_var <- length( vs_prob )
  full_samp <- n_var > 2
  
  betares <- matrix(0, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  vsres <- matrix(FALSE, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  colnames(betares) <- colnames(vsres) <- paste0( "beta_", 0L:(ncol(vsres)-1L) )
  
  ### First iteration
  betares[1L, ] <- start
  vsres[1L, ] <- varsel <- c( TRUE, logical(ncol(vsres) - 1) )
  
  X_mat <- X_mat_full[, varsel, drop = FALSE]
  V0i <- V0i_full[varsel, varsel]
  V0ib0 <- V0ib0_full[varsel]
  ldet_V0i <- sum( log(diag(chol(V0i))) )
  varsel_lprior <- log( dbinom(varsel, 1, vs_prob) )
  
  ### Main Loop
  for ( i in 2L:nrow(betares) ) {
    
    ### Propose dimension change
    change_ind <- ifelse( full_samp, sample(2L:n_var, 1L), 2L )    
    varsel_star <- vsres[i-1L, ]
    varsel_star[change_ind] <- !varsel_star[change_ind]
    
    X_mat_star <- X_mat_full[, varsel_star, drop = FALSE]
    V0i_star <- V0i_full[varsel_star, varsel_star]
    V0ib0_star <- V0ib0_full[varsel_star]
    V <- omega * nu
    pn <- psi * nu
    XtVi <- t( X_mat / V )
    XtVi_star <- t( X_mat_star / V )
    
    V_posti <- V0i + XtVi%*%X_mat
    V_posti_star <- V0i_star + XtVi_star%*%X_mat_star
    V_post <- chol2inv( chol(V_posti) )
    V_post_star <- chol2inv( chol(V_posti_star) )
    B_post <- V_post%*%( V0ib0 + XtVi%*%(z - pn) )
    B_post_star <- V_post_star%*%( V0ib0_star + XtVi_star%*%(z - pn) )
    
    ldet_V_post <- sum( log(diag(chol(V_post))) )
    ldet_V_post_star <- sum( log(diag(chol(V_post_star))) )
    ldet_V0i_star <- sum( log(diag(chol(V0i_star))) )
    varsel_lprior_star <- log( dbinom(varsel_star, 1, vs_prob) )
    
    lkernel <- 0.5*t(B_post)%*%V_posti%*%B_post
    lkernel_star <- 0.5*t(B_post_star)%*%V_posti_star%*%B_post_star
    ldenom <- -sum( ldet_V_post, ldet_V0i, lkernel, varsel_lprior )
    lnum <- sum( ldet_V_post_star, ldet_V0i_star, lkernel_star, varsel_lprior_star )
    
    lu <- log(runif(1))
    if ( lu <= sum(lnum, ldenom) ) {
      vsres[i, ] <- varsel_star
      X_mat <- X_mat_star
      V0i <- V0i_star
      V0ib0 <- V0ib0_star
      V_post <- V_post_star
      B_post <- B_post_star
    }
    else
      vsres[i, ] <- vsres[i-1L, ]
    
    ### Update beta
    beta <- betares[i, vsres[i, ]] <- mvrnorm( mu = B_post, Sigma = V_post )
    
    ### Update z
    Xb <- X_mat %*% beta
    z <- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
    
    ### Update nu
    mu <- sqrt( lambda * omega / (z - Xb)^2  )
    nu <- 1 / rinvgauss( n = n_y, mean = mu, shape = lambda )
    
    ### Print progress
    if ( !i%%500L )
      print( paste0("Current iteration: ", i) )
    
  }
  
  # Return results
  list( betadraw = betares, vardraw = vsres, MCMC_start_val = start,
        prior = prior, quantile = quantile, MCMC_length = MCMC_length,
        response_vec = yv, design_mat = X_mat_full )
  
}


