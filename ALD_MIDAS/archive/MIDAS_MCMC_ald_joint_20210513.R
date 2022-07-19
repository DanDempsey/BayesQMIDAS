###### MCMC Routines for Quantile IRTS MIDAS Inference 
###### Daniel Dempsey

# To be 
ins <- function(a, to.insert, pos) {

  if (pos <= 1)
    return( c(to.insert, a[seq(pos[1], length(a))]) )
  if (pos >= length(a)+1)
    return( c(a, to.insert)  )
  
  c(a[seq(pos[1]-1)], to.insert, a[seq(pos[1], length(a))])
  
}

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
      post_old <- (-0.5*(crossprod(Xb_old) - z_pn%*%Xb_old) + gam[j]) +
        dnorm(p_old[1], sd = 10, log = TRUE) + log(dtruncnorm(p_old[2], b=0, sd = 10))
      
      # Proposed posterior
      p_star <- DLF_star[pinds]
      Zenv$regime_object[[j]]$DLF_parameters <- p_star
      X_mat_star <- model.matrix( formula, Zenv )
      Xb_star <- X_mat_star %*% beta_new
      post_star <- (-0.5*(crossprod(Xb_star) - z_pn%*%Xb_star) + gam_star[j]) +
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
        response_vec = as.numeric(rtrunc), design_mat = X_mat )
  
}

MIDAS_MCMC_ald_vs <- function(formula, start, prior, quantile, 
                              adapt, MCMC_length, nburn, 
                              varsel_prior) {
  
  ### Prepare parameters
  Zenv <- environment( formula )
  formula2 <- formula( paste(c(as.character(formula)[2:1], 
                            "midas_design_matrices(regime_object2)"),
                            collapse = " "),
                          env = Zenv )
  formula3 <- formula( paste(c(as.character(formula)[2:1], 
                               "1"),
                             collapse = " "),
                       env = Zenv )
  rtrunc <- ifelse( Zenv$response_vector, TRUE, FALSE )
  n_y <- length(rtrunc)
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  lambda <- 2 + psi^2 / omega
  V0i_full <- chol2inv( chol(prior$V0) )
  V0i <- V0i_full[1, 1]
  V0ib0_full <- V0i_full %*% prior$beta0
  V0ib0 <- V0ib0_full[1, 1]
  
  X_mat_full <- model.matrix( formula, Zenv )
  X_mat <- model.matrix( formula3, Zenv )
  nu <- rep(1, n_y)
  V <- omega * nu
  
  # Beta initialization
  betares <- matrix(0, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  vsres <- matrix(FALSE, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  colnames(betares) <- colnames(vsres) <- paste0( "beta_", 0L:(ncol(vsres)-1L) )
  betares[1L, ] <- start
  beta_new <- start[1]
  new_name <- c()
  vsres[, 1] <- TRUE
  varsel <- varsel_star <- vsres[1, ]
  Xb <- X_mat %*% beta_new
  
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
  split_DLF <- split( 1:nDLF_sum, rep( 1:nvar, nDLF ) )
  
  accept <- rep( 0L, nvar )
  runs <- rep( 1L, nvar )
  sigmaDLF <- rep( 0.2, nvar )
  names(sigmaDLF) <- names(runs) <- names(accept) <- names(split_DLF) <- names(Zenv$regime_object)
  
  # Initialise covariate indicator
  if (missing(varsel_prior)) {
    varsel_prior <- c(1, rep(1/nvar, nvar))
  }
  varsel_lprior <- sum( dbinom(varsel, 1, varsel_prior, log = TRUE) )
  
  ### Main Loop
  for ( i in 2L:nrow(betares) ) {
    
    #if (i == 4) { browser() }
    
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
    V_post <- chol2inv( chol(V0i_full[varsel, varsel] + XtVi%*%X_mat) )
    B_post <- V_post%*%( V0ib0 + XtVi%*%z_pn/2 )
    beta_new <- mvrnorm( mu = B_post, Sigma = V_post )
    betares[i, varsel] <- beta_new
    names(beta_new) <- c("Int", new_name)
    
    ### Update DLF parameters
    
    # Pick out selected variables
    Zenv$regime_object2 <- Zenv$regime_object[varsel[-1]]
    
    # Variance adaptation
    for (j in seq_along(Zenv$regime_object2)) {
      
      #if (i == 10) { browser() }
      
      nn <- names(Zenv$regime_object2[j])
      runs_j <- runs[nn]
      accept_j <- accept[nn]
      
      if( adapt & (i < nburn) & (runs_j%%10 == 0) ){
        Delta <- 1/sqrt(runs_j)
        adapt_indicator <- ifelse( accept_j/runs_j > .234, TRUE, FALSE ) 
        if (adapt_indicator)
          sigmaDLF[nn] <- sigmaDLF[nn] * exp( Delta )
        else
          sigmaDLF[nn] <- sigmaDLF[nn] * exp( -Delta )
      }
      
      sigmaDLF_j <- sigmaDLF[nn]
      DLF_last <- Zenv$regime_object2[[j]]$DLF_parameters
      DLF_last[2] <- log(-DLF_last[2])
      runs[nn] <- runs[nn] + 1
      
      # Proposal
      DLF_star_gam <- rnorm(2, DLF_last, sigmaDLF_j)
      gam_star <- DLF_star_gam[2]
      DLF_star <- DLF_star_gam
      DLF_star[2] <- -exp(DLF_star[2])
      
      pinds <- split_DLF[[nn]]
      
      # Current posterior
      p_old <- Zenv$regime_object2[[j]]$DLF_parameters
      Xb_old <- X_mat %*% beta_new
      post_old <- (-0.5*(crossprod(Xb_old) - z_pn%*%Xb_old) + DLF_last[2]) +
        dnorm(DLF_last[1], sd = 10, log = TRUE) + 
        log(dtruncnorm(-exp(DLF_last[2]), b = 0, sd = 10))
      
      # Proposed posterior
      #p_star <- DLF_star[pinds]
      Zenv$regime_object2[[j]]$DLF_parameters <- DLF_star
      X_mat_star <- model.matrix( formula2, Zenv )
      Xb_star <- X_mat_star %*% beta_new
      post_star <- (-0.5*(crossprod(Xb_star) - z_pn%*%Xb_star) + gam_star) +
        dnorm(DLF_star[1], sd = 10, log = TRUE) + 
        log(dtruncnorm(DLF_star[2], b = 0, sd = 10))
      
      acp <- post_star - post_old
      if (acp >= log(runif(1))) {
        DLFres[[nn]][i, ] <- DLF_star
        gam <- gam_star
        X_mat <- X_mat_star
        Xb <- Xb_star
        accept[j] <- accept[j] + 1L
        DLFres_mat[i, pinds] <- DLF_star_gam
      }
      else {
        DLFres[[j]][i, ] <- p_old
        Zenv$regime_object2[[j]]$DLF_parameters <- p_old
        DLFres_mat[i, pinds] <- DLFres_mat[i-1, pinds]
      }
    }
    Zenv$regime_object[varsel[-1]] <- Zenv$regime_object2[varsel[-1]]
    
    ### Update Covariance Indicator
    if (i == 7) { browser() }
    
    lDLF_prior <- 0
    for (k in seq_along(Zenv$regime_object2)) {
      DLF_pars <- Zenv$regime_object2[[k]]$DLF_parameters
      lDLF_prior <- sum(lDLF_prior, dnorm(DLF_pars[1], sd = 10, log = TRUE), 
                        log(dtruncnorm(DLF_pars[2], b = 0, sd = 10)))
    }
    
    lpost <- sum( -0.5*(crossprod(Xb) - z_pn%*%Xb), 
                 dnorm(beta_new, 0, 10, log = TRUE),
                 sum(dbinom(varsel, 1, varsel_prior, log = TRUE)),
                 lDLF_prior )
    change_ind <- sample(1L:nvar, 1L) + 1L
    varsel_star[change_ind] <- !varsel_star[change_ind]
    new_name2 <- names(Zenv$regime_object)[varsel_star[1:nvar + 1]]
    change_name <- names(Zenv$regime_object)[change_ind-1]
    
    if( varsel[change_ind] ) {
      # Move to a lower dimension
      name_drop <- names(Zenv$regime_object[change_ind-1L])
      Theta_drop <- c(beta_new[name_drop], 
                      Zenv$regime_object2[[name_drop]]$DLF_parameters)
      T_drop <- Theta_drop
      T_drop <- log(-T_drop[3])
      beta_drop <- beta_new[setdiff(names(beta_new), name_drop)]
      
      Zenv$regime_object2 <- Zenv$regime_object[new_name2]
      X_mat_star <- model.matrix( formula2, Zenv )
      Xb_star <- X_mat_star %*% beta_drop
      
      lDLF_prior_star <- 0
      for (k in seq_along(Zenv$regime_object2)) {
        DLF_pars <- Zenv$regime_object2[[k]]$DLF_parameters
        lDLF_prior_star <- sum(lDLF_prior_star, dnorm(DLF_pars[1], sd = 10, log = TRUE), 
                               log(dtruncnorm(DLF_pars[2], b=0, sd = 10)))
      }
      
      lpost_star <- sum( -0.5*(crossprod(Xb_star) - z_pn%*%Xb_star),
                        dnorm(beta_drop, 0, 10, log = TRUE),
                        sum(dbinom(varsel_star, 1, varsel_prior, log = TRUE)),
                        lDLF_prior_star )
      
      lproposal <- sum(dnorm(T_drop, 0, 0.2, log = TRUE))
      ljacob <- abs( log(-Theta_drop[3]) )
      #accept_prob <- sum( lpost_star, -lpost,
      #                    -lproposal, ljacob )
      accept_prob <- sum( lpost_star, -lpost, lproposal, ljacob )
      
      if (accept_prob >= log( runif(1) )) {
        #browser()
        varsel <- vsres[i, ] <- varsel_star
        beta_new <- beta_drop
        X_mat <- X_mat_star
        Xb <- Xb_star
        Zenv$regime_object[[change_ind-1]]$DLF_parameters <- NULL
        new_name <- new_name2
      }
      else {
        vsres[i, ] <- vsres[i - 1L, ]
        varsel_star <- varsel
      }
    }
    else {
      # Move to a higher dimension
      Zenv$regime_object2 <- Zenv$regime_object[new_name2]
      beta_star_new <- rnorm(1, 0, 0.2)
      beta_star <- ins(beta_new, beta_star_new, change_ind)
      names(beta_star) <- c("Int", new_name2)
      DLF_star <- rnorm(2, 0, 0.2)
      gam_star <- DLF_star[2]
      DLF_star[2] <- -exp(DLF_star[2])
      
      Zenv$regime_object2[[change_name]]$DLF_parameters <- DLF_star
      X_mat_star <- model.matrix( formula2, Zenv )
      Xb_star <- X_mat_star %*% beta_star 
      
      lDLF_prior_star <- 0
      for (k in seq_along(Zenv$regime_object2)) {
        DLF_pars <- Zenv$regime_object2[[k]]$DLF_parameters
        lDLF_prior_star <- sum(lDLF_prior_star, dnorm(DLF_pars[1], sd = 10, log = TRUE), 
                               log(dtruncnorm(DLF_pars[2], b=0, sd = 10)))
      }
      
      lpost_star <- sum( -0.5*(crossprod(Xb_star) - z_pn%*%Xb_star),
                        dnorm(beta_star, 0, 10, log = TRUE),
                        sum(dbinom(varsel_star, 1, varsel_prior, log = TRUE)),
                        lDLF_prior_star )
      
      lproposal <- sum(dnorm(c(beta_star_new, DLF_star), 0, 0.2, log = TRUE))
      ljacob <- gam_star
      accept_prob <- sum( lpost_star, -lpost,
                          -lproposal, ljacob )
      
      if (accept_prob >= log( runif(1) )) {
        varsel <- vsres[i, ] <- varsel_star
        beta_new <- beta_star
        X_mat <- X_mat_star
        Xb <- Xb_star
        Zenv$regime_object[[change_ind-1]]$DLF_parameters <- DLF_star
        new_name <- new_name2
      }
      else {
        vsres[i, ] <- vsres[i - 1L, ]
        varsel_star <- varsel
      }
    }
    
    ### Print progress
    if ( !i%%500L )
      print( paste0("Current iteration: ", i) )

  }
  
  list( betadraw = betares, DLFres = DLFres, varsel = vsres, 
        MCMC_start_val = start, prior = prior, quantile = quantile, 
        MCMC_length = MCMC_length, response_vec = as.numeric(rtrunc), 
        design_mat = X_mat )
  
}


