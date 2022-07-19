##### MCMC Component Functions
##### Daniel Dempsey

# Helper function to insert a value inside a vector
#ins <- function(a, to.insert, pos) {
#
#  pos <- pos + 1L
#  if (pos <= 1L)
#    return( c(to.insert, a[seq(pos[1L], length(a))]) )
#  if (pos >= length(a)+1L)
#    return( c(a, to.insert)  )
#  
#  c(a[seq(pos[1L]-1L)], to.insert, a[seq(pos[1L], length(a))])
#  
#}
#
## Function to build design matrix
#make_design_mat <- function(regime_object, offset_only, vars_ind) {
#  
#  if ( !any(vars_ind) ) 
#    return( offset_only )
#  
#  cbind( 1, midas_design_matrices(regime_object[vars_ind]) )
#  
#}

#all_DLF_prior <- function(x, m = c(0, 0), s = c(1, 1)) {
#  sum( dnorm(x[1L], mean = m[1], sd = s[1], log = TRUE), 
#       log(dtruncnorm(x[2L], b = 0, mean = m[2], sd = s[2])) )
#}

MCMC_setup_vs <- function(envir_object, start, prior, nburn, adapt,
                          varsel_prior, quantile, MCMC_length) {
  
  regime_object <- envir_object$regime_object
  rtrunc <- ifelse( envir_object$response_vector, TRUE, FALSE )
  n_y <- length(rtrunc)
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  lambda <- 2 + psi^2 / omega
  V0i_full <- chol2inv( chol(prior$V0) )
  V0i <- V0i_full[1L, 1L]
  V0ib0_full <- V0i_full %*% prior$beta0
  V0ib0 <- V0ib0_full[1L, 1L]
  
  int <- matrix(1, nrow = n_y, ncol = 1)
  varsel <- rep( FALSE, length(regime_object) )
  X_mat_full <- make_design_mat( regime_object, int, !varsel )
  X_mat <- make_design_mat( regime_object, int, varsel )
  nu <- rep(1, n_y)
  V <- omega * nu
  z <- rep(0, n_y)
  pn <- psi * nu
  z_pn <- z - pn
  
  # Beta initialization
  betares <- matrix(0, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  vsres <- matrix(FALSE, ncol = ncol(X_mat_full)-1L, nrow = MCMC_length + 1L)
  varnames <<- names(envir_object$regime_object)
  colnames(betares) <- c( "Int", varnames )
  colnames(vsres) <- varnames
  betares[1L, ] <- start
  varsel <- vsres[1L, ]
  Xb <- X_mat %*% start[1L]
  
  # DLF parameter initialization
  var_extract <- paste0('regime_object$', varnames, '$DLF_parameters')
  DLFres <- Map( parse, text = var_extract )
  DLFres <- lapply(DLFres, eval, envir = envir_object)
  names(DLFres) <- varnames
  make_mat <- function(x) {
    matrix(rep(x, each = MCMC_length + 1L), nrow = MCMC_length + 1L)
  }
  DLFres <- lapply(DLFres, make_mat)
  
  nvar <- ncol( betares ) - 1L
  varvec <- 1L:nvar
  accept <- accept_start <- runs <- runs_start <- rep( 1L, nvar )
  sigmaDLF <- sigmaDLF_start <- rep( 0.2, nvar )
  names(sigmaDLF) <- names(runs) <- names(accept) <- names(runs_start) <- names(sigmaDLF_start) <- varnames
  
  #browser()
  ### Return variables
  rm(var_extract, make_mat, X_mat_full, nvar, start, z)
  environment()
  
}

update_covar <- function(i) {
  
  ### Current
  varsel_int <- c( TRUE, varsel )
  V0i <- V0i_full[varsel_int, varsel_int]
  V0ib0 <- V0ib0_full[varsel_int]
  ldet_V0i <- sum( log(diag(chol(V0i))) )
  varsel_lprior <- dbinom(varsel, 1, varsel_prior, log = TRUE)
  ldet_V0i <- sum( log(diag(chol(V0i))) )
  XtVi <- t( X_mat / V )
  
  V_posti <- V0i + XtVi%*%X_mat
  V_post <<- V_post <- chol2inv( chol(V_posti) )
  B_post <<- B_post <- V_post%*%( V0ib0 + XtVi%*%z_pn )
  
  if ( i%%10 != 0 ) { 
    vsres[i, ] <<- vsres[i-1, ]
    return() 
  }
  
  ### Propose dimension change
  change_ind <- sample( varvec, 1L )    
  varsel_star <- vsres[i-1L, ]
  varsel_star[change_ind] <- !varsel_star[change_ind]
  varsel_star_int <- c( TRUE, varsel_star )
  change_name <- names(regime_object)[change_ind]
  
  if( varsel[change_ind] ) {
    
    # Move to a lower dimension
    Theta_drop <- regime_object[[change_name]]$DLF_parameters
    lproposal <- dnorm(c(Theta_drop[1], log(-Theta_drop[2])), 0, sigmaDLF_start[change_name], log = TRUE)
    ljacob <- -2 * log( -Theta_drop[2] )
    lDLFprior <- all_DLF_prior( Theta_drop )
    DLF_component <- sum( lproposal, ljacob, -lDLFprior )
    
  }
  else {
    
    # Proposal variance adaptation
    runs_j <- runs_start[change_name]
    #if( adapt & (i < nburn) & (runs_j%%200L == 0L) ){
    #  Delta <- 1/sqrt(runs_j)
    #  adapt_indicator <- ifelse( accept_start[change_name]/runs_j > .234, TRUE, FALSE ) 
    #  if (adapt_indicator)
    #    sigmaDLF_start[change_name] <<- sigmaDLF_start[change_name] * exp( Delta )
    #  else
    #    sigmaDLF_start[change_name] <<- sigmaDLF_start[change_name] * exp( -Delta )
    #}
    
    #sigmaDLF_j <- sigmaDLF_start[change_name]
    sigmaDLF_j <- 0.2
    runs_start[change_name] <<- runs_start[change_name] + 1L
    
    # Move to a higher dimension
    Theta_star <- rnorm(2L, 0, sigmaDLF_j)
    lproposal <- dnorm(Theta_star, 0, sigmaDLF_j, log = TRUE)
    ljacob <- 2*Theta_star[2]
    Theta_star[2] <- -exp(Theta_star[2])
    regime_object[[change_name]]$DLF_parameters <- Theta_star
    lDLFprior <- all_DLF_prior( Theta_star )
    DLF_component <- sum( -lproposal, ljacob, lDLFprior )
    
  }
  
  X_mat_star <- make_design_mat( regime_object, int, varsel_star )
  V0i_star <- V0i_full[varsel_star_int, varsel_star_int]
  V0ib0_star <- V0ib0_full[varsel_star_int]
  XtVi_star <- t( X_mat_star / V )
  
  V_posti_star <- V0i_star + XtVi_star%*%X_mat_star
  V_post_star <- chol2inv( chol(V_posti_star) )
  B_post_star <- V_post_star%*%( V0ib0_star + XtVi_star%*%z_pn )
  
  ldet_V_post <- sum( log(diag(chol(V_post))) )
  ldet_V_post_star <- sum( log(diag(chol(V_post_star))) )
  ldet_V0i_star <- sum( log(diag(chol(V0i_star))) )
  varsel_lprior_star <- dbinom(varsel_star, 1, varsel_prior, log = TRUE)
  
  lkernel <- 0.5*crossprod(B_post, V_posti)%*%B_post
  lkernel_star <- 0.5*crossprod(B_post_star, V_posti_star)%*%B_post_star
  ldenom <- -sum( ldet_V_post, ldet_V0i, lkernel, varsel_lprior )
  lnum <- sum( ldet_V_post_star, ldet_V0i_star, lkernel_star, 
               varsel_lprior_star, DLF_component )
  
  #print(sum(lnum, ldenom))
  if ( log(runif(1L)) <= sum(lnum, ldenom) ) {
    if( !varsel[change_ind] ) {
      regime_object[[change_name]]$DLF_parameters <<- Theta_star
      accept_start[change_name] <<- accept_start[change_name] + 1L
    }
    vsres[i, ] <<- varsel <<- varsel_star
    X_mat <<- X_mat_star
    V0i <<- V0i_star
    V0ib0 <<- V0ib0_star
    V_post <<- V_post_star
    B_post <<- B_post_star
  }
  else
    vsres[i, ] <<- vsres[i-1L, ]
  
}

update_qr <- function(i) {
  
  ### Update beta
  bet <- betares[i, c(TRUE, vsres[i, ])] <<- mvrnorm( mu = B_post, Sigma = V_post )
  
  ### Update z
  Xb <<- Xb <- X_mat %*% bet
  z <- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
  
  ### Update nu
  mu <- sqrt( lambda * omega / (z - Xb)^2  )
  nu <- 1 / rinvgauss( n = n_y, mean = mu, shape = lambda )
  z_pn <<- z - psi*nu
  Xb <<- Xb
  
}

update_dlf <- function(i) {
  
  ll <- -0.5*crossprod(Xb) + z_pn%*%Xb
  bet <- betares[i, c(TRUE, varsel)]
  
  # Loop over DLF updates
  for (j in which(varsel)) {
    
    #browser()
    # Proposal variance adaptation
    if( adapt & (i < nburn) & (runs[j]%%200L == 0L) ){
      Delta <- 1/sqrt(runs[j])
      adapt_indicator <- ifelse( accept[j]/runs[j] > .234, TRUE, FALSE ) 
      if (adapt_indicator)
        sigmaDLF[j] <<- sigmaDLF[j] * exp( Delta )
      else
        sigmaDLF[j] <<- sigmaDLF[j] * exp( -Delta )
    }
    runs[j] <<- runs[j] + 1L
    
    # Current posterior
    DLF_last <- DLF_trans_last <- regime_object[[j]]$DLF_parameters
    post_old <- ll + all_DLF_prior(DLF_last)
    
    # Proposal
    DLF_trans_last[2L] <- log( -DLF_trans_last[2L] )
    DLF_star <- rnorm( 2L, DLF_trans_last, sigmaDLF[j] )
    ljacob_star <- 2 * ( DLF_star[2L] - DLF_trans_last[2L] )
    DLF_star[2L] <- -exp( DLF_star[2L] )
    
    # Proposed posterior
    regime_object[[j]]$DLF_parameters <- DLF_star
    Xb_star <- cbind( 1, midas_design_matrices(regime_object[varsel]) ) %*% bet
    ll_star <- -0.5*crossprod(Xb_star) + z_pn%*%Xb_star
    post_star <- ll_star + all_DLF_prior(DLF_star) + ljacob_star
    
    if ( (post_star - post_old) >= log(runif(1)) ) {
      DLFres[[j]][i, ] <<- regime_object[[j]]$DLF_parameters <<- DLF_star
      accept[j] <<- accept[j] + 1L
      ll <- ll_star
      Xb <<- Xb_star
    }
    else {
      DLFres[[j]][i, ] <<- DLF_last
    }
  }
  
}

