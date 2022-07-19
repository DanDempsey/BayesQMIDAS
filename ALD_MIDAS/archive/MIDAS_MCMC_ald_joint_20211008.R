###### MCMC Routines for Quantile IRTS MIDAS Inference 
###### Daniel Dempsey

MIDAS_MCMC_ald_fixedDLF <- function(formula, start, prior, quantile, MCMC_length) {
  
  ### Prepare parameters
  Zenv <- environment( formula )
  rtrunc <- ifelse( Zenv$response_vector, TRUE, FALSE )
  n_y <- length(rtrunc)
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  lambda <- 2 + psi^2 / omega
  V0i <- chol2inv( chol(prior$V0) )
  V0ib0 <- V0i %*% prior$beta0
  
  X_mat <- model.matrix( formula, Zenv )
  nu <- rep(1, n_y)
  pn <- psi * nu
  V <- omega * nu
  
  betares <- matrix(0, ncol = ncol(X_mat), nrow = MCMC_length + 1L)
  colnames(betares) <- paste0( "beta_", 0L:(ncol(betares)-1L) )
  betares[1L, ] <- start$beta

  ### Main Loop
  for ( i in 2L:nrow(betares) ) {
    
    ### Update z
    Xb <- X_mat %*% betares[i-1L, ]
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
    betares[i, ] <- mvrnorm( mu = B_post, Sigma = V_post )
    
    ### Print progress
    if ( !i%%500L )
      print( paste0("Current iteration: ", i) )

  }
  
  list( betadraw = betares, MCMC_start_val = start, prior = prior, 
        quantile = quantile, MCMC_length = MCMC_length,
        response_vec = as.numeric(rtrunc), design_mat = X_mat_full )
  
}

MIDAS_MCMC_ald_fixedDLF_vs <- function(formula, start, prior, quantile, MCMC_length) {
  
  ### Preparation
  Zenv <- environment( formula )
  rtrunc <- ifelse( Zenv$response_vector, TRUE, FALSE )
  n_y <- length(rtrunc)
  
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
  varvec <- seq(2, n_var)
  
  betares <- matrix(0, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  vsres <- matrix(FALSE, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  colnames(betares) <- colnames(vsres) <- paste0( "beta_", 0L:(ncol(vsres)-1L) )
  
  ### First iteration
  betares[1L, ] <- start$beta
  vsres[1L, ] <- varsel <- c( TRUE, logical(ncol(vsres) - 1) )
  
  X_mat <- X_mat_full[, varsel, drop = FALSE]
  V0i <- V0i_full[varsel, varsel]
  V0ib0 <- V0ib0_full[varsel]
  ldet_V0i <- sum( log(diag(chol(V0i))) )
  varsel_lprior <- log( dbinom(varsel, 1, vs_prob) )
  
  ### Main Loop
  for ( i in 2L:nrow(betares) ) {
    
    ### Propose dimension change
    change_ind <- sample( varvec, 1L )    
    varsel_star <- vsres[i-1L, ]
    varsel_star[change_ind] <- !varsel_star[change_ind]
    
    X_mat_star <- X_mat_full[, varsel_star, drop = FALSE]
    V0i_star <- V0i_full[varsel_star, varsel_star]
    V0ib0_star <- V0ib0_full[varsel_star]
    V <- omega * nu
    pn <- psi * nu
    z_pn <- z - pn
    XtVi <- t( X_mat / V )
    XtVi_star <- t( X_mat_star / V )
    
    V_posti <- V0i + XtVi%*%X_mat
    V_posti_star <- V0i_star + XtVi_star%*%X_mat_star
    V_post <- chol2inv( chol(V_posti) )
    V_post_star <- chol2inv( chol(V_posti_star) )
    B_post <- V_post%*%( V0ib0 + XtVi%*%z_pn )
    B_post_star <- V_post_star%*%( V0ib0_star + XtVi_star%*%z_pn )
    
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
        response_vec = as.numeric(rtrunc), design_mat = X_mat_full )
  
}

MIDAS_MCMC_ald <- function(formula, start, prior, quantile, 
                           adapt, MCMC_length, nburn) {
  
  ### Prepare parameters
  Zenv <- environment( formula )
  rtrunc <- ifelse( Zenv$response_vector, TRUE, FALSE )
  n_y <- length(rtrunc)
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  lambda <- 2 + psi^2 / omega
  V0i <- chol2inv( chol(prior$V0) )
  V0ib0 <- V0i %*% prior$beta0
  
  X_mat <- model.matrix( formula, Zenv )
  nu <- rep(1, n_y)
  V <- omega * nu
  
  # Beta initialization
  betares <- matrix(0, ncol = ncol(X_mat), nrow = MCMC_length + 1L)
  colnames(betares) <- paste0( "beta_", 0L:(ncol(betares)-1L) )
  betares[1L, ] <- beta_new <- start
  
  # DLF parameter initialization
  varnames <- names(Zenv$regime_object)
  var_extract <- paste0('regime_object$', varnames, '$DLF_parameters')
  DLFres <- Map( parse, text = var_extract )
  DLFres <- lapply(DLFres, eval, envir = Zenv)
  names(DLFres) <- varnames
  make_mat <- function(x) {
    matrix(rep(x, each = MCMC_length + 1L), nrow = MCMC_length + 1L)
  }
  DLFres <- lapply(DLFres, make_mat)
  DLFres_mat <- do.call('cbind', DLFres)
  rm(varnames, var_extract, make_mat)
  
  nvar <- ncol( betares ) - 1
  nDLF <- sapply( DLFres, ncol )
  nDLF_sum <- sum( nDLF )
  
  trans_ind <- cumsum( nDLF )
  DLFres_mat[, trans_ind] <- log(-DLFres_mat[, trans_ind])
  gam <- DLFres_mat[1, trans_ind]
  split_DLF <- split(1:nDLF_sum, rep( 1:nvar, nDLF ) )
  
  accept <- rep( 0L, nvar )
  sigmaDLF <- rep( 0.2, nvar )
  
  ### Main Loop
  for ( i in 2L:nrow(betares) ) {
    
    ### Update z
    Xb <- X_mat %*% beta_new
    z <- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
    
    ### Update nu
    mu <- sqrt( lambda * omega / (z - Xb)^2  )
    nu <- 1 / rinvgauss( n = n_y, mean = mu, shape = lambda )
    
    ### Update beta
    V <- omega * nu
    z_pn <- 2*(z - psi*nu)
    XtVi <- t( X_mat / V )
    V_post <- chol2inv( chol(V0i + XtVi%*%X_mat) )
    B_post <- V_post%*%( V0ib0 + XtVi%*%z_pn/2 )
    beta_new <- mvrnorm( mu = B_post, Sigma = V_post )
    betares[i, ] <- beta_new
    
    ### Update DLF parameters
    # Variance adaptation
    if( adapt & (i < nburn) & (i%%10 == 0) ){
      Delta <- 1/sqrt(i)
      #Delta <- ifelse( delta < .01, delta, .01 )
      adapt_ind <- ifelse( accept/(i-1) > .234, TRUE, FALSE ) 
      sigmaDLF[adapt_ind] <- sigmaDLF[adapt_ind] * exp( Delta ) 
      sigmaDLF[!adapt_ind] <- sigmaDLF[!adapt_ind] * exp( -Delta )
    }
    
    # Proposal
    DLF_star_gam <- rnorm(nDLF_sum, DLFres_mat[i-1, ], rep(sigmaDLF, each = 2L))
    gam_star <- DLF_star_gam[trans_ind]
    DLF_star <- DLF_star_gam
    DLF_star[trans_ind] <- -exp(DLF_star[trans_ind])
    
    for (j in 1:nvar) {
      pinds <- split_DLF[[j]]
      
      # Current posterior
      p_old <- Zenv$regime_object[[j]]$DLF_parameters
      Xb_old <- X_mat %*% beta_new
      post_old <- -0.5*(crossprod(Xb_old) - z_pn%*%Xb_old) + #gam[j] +
        dnorm(p_old[1], sd = 10, log = TRUE) + log(dtruncnorm(p_old[2], b=0, sd = 10))
      
      # Proposed posterior
      p_star <- DLF_star[pinds]
      Zenv$regime_object[[j]]$DLF_parameters <- p_star
      X_mat_star <- model.matrix( formula, Zenv )
      Xb_star <- X_mat_star %*% beta_new
      post_star <- -0.5*(crossprod(Xb_star) - z_pn%*%Xb_star) + gam_star[j] - log( -p_old[2L] ) +
        dnorm(p_star[1], sd = 10, log = TRUE) + log(dtruncnorm(p_star[2], b=0, sd = 10))
      
      acp <- min( 0, post_star - post_old )
      if (acp >= log(runif(1))) {
        DLFres[[j]][i, ] <- p_star
        gam <- gam_star
        X_mat <- X_mat_star
        accept[j] <- accept[j] + 1L
        DLFres_mat[i, pinds] <- DLF_star_gam[pinds]
      }
      else {
        DLFres[[j]][i, ] <- p_old
        Zenv$regime_object[[j]]$DLF_parameters <- p_old
        DLFres_mat[i, pinds] <- DLFres_mat[i-1, pinds]
      }
    }
    
    ### Print progress
    if ( !i%%500L )
      print( paste0("Current iteration: ", i) )

  }
  
  list( betadraw = betares, DLFres = DLFres, MCMC_start_val = start, 
        prior = prior, quantile = quantile, MCMC_length = MCMC_length,
        response_vec = as.numeric(rtrunc), design_mat = X_mat,
        acceptance = accept )
  
}

#MIDAS_MCMC_ald_vs <- function(formula, start, prior, quantile, 
#                              adapt, MCMC_length, nburn, 
#                              varsel_prior) {
#  
#  ### Prepare data environment
#  Zenv <- environment( formula )
#  if (missing(varsel_prior)) {
#    nvar <- length(start) - 1
#    varsel_prior <- rep(1/nvar, nvar)
#  }
#  dat_ee <- MCMC_setup_vs(envir_object = Zenv, start = start, prior = prior, 
#                          varsel_prior = varsel_prior, quantile = quantile, 
#                          MCMC_length = MCMC_length, nburn = nburn, 
#                          adapt = adapt)
#  parent.env(dat_ee) <- globalenv()
#  
#  ### Set environments
#  environment(update_covar) <- environment(update_qr) <- environment(update_dlf) <- dat_ee
#  
#  ### Main Loop
#  for ( i in 2L:(MCMC_length+1L) ) {
#    
#    update_covar(i)
#    update_qr(i)
#    update_dlf(i)
#    
#    # Print progress
#    if ( !i%%500L )
#      print( paste0("Current iteration: ", i) )
#
#  }
#  
#  list(betares = dat_ee$betares, DLFres = dat_ee$DLFres, vsres = dat_ee$vsres,
#       accept = dat_ee$accept, runs = dat_ee$runs, sigmaDLF = dat_ee$sigmaDLF,
#       accept_start = dat_ee$accept_start, runs_start = dat_ee$runs_start, 
#       sigmaDLF_start = dat_ee$sigmaDLF_start)
#  
#}
#
#
#