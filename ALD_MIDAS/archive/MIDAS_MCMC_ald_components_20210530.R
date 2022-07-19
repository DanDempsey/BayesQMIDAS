##### MCMC Component Functions
##### Daniel Dempsey

# Helper function to insert a value inside a vector
ins <- function(a, to.insert, pos) {

  pos <- pos + 1L
  if (pos <= 1L)
    return( c(to.insert, a[seq(pos[1L], length(a))]) )
  if (pos >= length(a)+1L)
    return( c(a, to.insert)  )
  
  c(a[seq(pos[1L]-1L)], to.insert, a[seq(pos[1L], length(a))])
  
}

# Function to build design matrix
make_design_mat <- function(regime_object, n_y, vars_ind) {
  
  if ( !any(vars_ind) ) 
    return( matrix(1, nrow = n_y, ncol = 1) )
  
  cbind( 1, midas_design_matrices(regime_object[vars_ind]) )
  
}

all_DLF_prior <- function(k) {
  sum(dnorm(k$DLF_parameters[1L], sd = 10, log = TRUE), 
      log(dtruncnorm(k$DLF_parameters[2L], b = 0, sd = 10)))
}

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
  
  varsel <- rep( FALSE, length(regime_object) )
  X_mat_full <- make_design_mat( regime_object, n_y, !varsel )
  X_mat <- make_design_mat( regime_object, n_y, varsel )
  nu <- rep(1, n_y)
  V <- omega * nu
  
  # Beta initialization
  betares <- matrix(0, ncol = ncol(X_mat_full), nrow = MCMC_length + 1L)
  vsres <- matrix(FALSE, ncol = ncol(X_mat_full)-1L, nrow = MCMC_length + 1L)
  varnames <<- names(envir_object$regime_object)
  colnames(betares) <- c( "Int", varnames )
  colnames(vsres) <- varnames
  betares[1L, ] <- start
  beta_new <- start[1L]
  varsel <- vsres[1L, ]
  Xb <- X_mat %*% beta_new
  
  # DLF parameter initialization
  var_extract <- paste0('regime_object$', varnames, '$DLF_parameters')
  DLFres <- Map( parse, text = var_extract )
  DLFres <- lapply(DLFres, eval, envir = envir_object)
  names(DLFres) <- varnames
  make_mat <- function(x) {
    matrix(rep(x, each = MCMC_length + 1L), nrow = MCMC_length + 1L)
  }
  DLFres <- lapply(DLFres, make_mat)
  rm(var_extract, make_mat)
  
  nvar <- ncol( betares ) - 1
  accept <- rep( 0L, nvar )
  runs <- rep( 1L, nvar )
  sigmaDLF <- rep( 0.2, nvar )
  names(sigmaDLF) <- names(runs) <- names(accept) <- varnames
  
  ### Return variables
  environment()
  
}

# Quantile Regression Section

update_qr <- function(i) {
  
  ### Update z
  Xb <<- X_mat %*% beta_new
  z <<- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
  
  ### Update nu
  mu <<- sqrt( lambda * omega / (z - Xb)^2  )
  nu <<- 1 / rinvgauss( n = n_y, mean = mu, shape = lambda )
  
  ### Update beta
  varsel_int <- c( TRUE, varsel )
  V <<- omega * nu
  z_pn <<- 2*(z - psi*nu)
  XtVi <<- t( X_mat / V )
  V_post <<- chol2inv( chol(V0i_full[varsel_int, varsel_int] + XtVi%*%X_mat) )
  B_post <<- V_post%*%( V0ib0 + XtVi%*%z_pn/2 )
  beta_new <<- betares[i, varsel_int] <<- mvrnorm( mu = B_post, Sigma = V_post )
  names(beta_new) <<- c( "Int", varnames[varsel] )
  
}

update_dlf <- function(i) {
  
  ll <<- -0.5*(crossprod(Xb) - z_pn%*%Xb)
  
  # Loop over DLF updates
  for (j in which(varsel)) {
    
    runs_j <- runs[j]
    accept_j <- accept[j]
    
    # Variance adaptation
    if( adapt & (i < nburn) & (runs_j%%10L == 0L) ){
      Delta <- 1/sqrt(runs_j)
      adapt_indicator <- ifelse( accept_j/runs_j > .234, TRUE, FALSE ) 
      if (adapt_indicator)
        sigmaDLF[j] <<- sigmaDLF[j] * exp( Delta )
      else
        sigmaDLF[j] <<- sigmaDLF[j] * exp( -Delta )
    }
    
    sigmaDLF_j <- sigmaDLF[j]
    runs[j] <<- runs[j] + 1L
    
    # Current posterior
    DLF_last <- p_old <- regime_object[[j]]$DLF_parameters
    DLF_last[2L] <- log(-DLF_last[2L])
    post_old <- ll +
      dnorm(DLF_last[1L], sd = 10, log = TRUE) + 
      log(dtruncnorm(p_old[2L], b = 0, sd = 10)) +
      p_old[2L]
    
    # Proposal
    DLF_star <- rnorm(2L, DLF_last, sigmaDLF_j)
    DLF_star[2L] <- -exp(DLF_star[2L])
    
    # Proposed posterior
    regime_object[[j]]$DLF_parameters <<- DLF_star
    X_mat_star <- make_design_mat( regime_object, n_y, varsel )
    Xb_star <- X_mat_star %*% beta_new
    ll_star <- -0.5*(crossprod(Xb_star) - z_pn%*%Xb_star)
    post_star <- ll_star +
      dnorm(DLF_star[1], sd = 10, log = TRUE) + 
      log(dtruncnorm(DLF_star[2], b = 0, sd = 10)) +
      DLF_star[2L]
    
    acp <- post_star - post_old
    if (acp >= log(runif(1))) {
      DLFres[[j]][i, ] <<- DLF_star
      X_mat <<- X_mat_star
      ll <<- ll_star
      accept[j] <<- accept[j] + 1L
    }
    else {
      regime_object[[j]]$DLF_parameters <<- DLFres[[j]][i, ] <<- p_old
    }
  }
  
}

update_covars <- function(i) {
  
  lDLF_prior <- 0
  for (k in regime_object[varsel]) {
    DLF_pars <- k$DLF_parameters
    lDLF_prior <- sum(lDLF_prior, dnorm(DLF_pars[1L], sd = 10, log = TRUE), 
                      log(dtruncnorm(DLF_pars[2L], b = 0, sd = 10)))
  }
  
  lpost <- sum( ll, 
               dnorm(beta_new, 0, 10, log = TRUE),
               sum(dbinom(varsel, 1, varsel_prior, log = TRUE)),
               lDLF_prior )
  change_ind <- sample(1L:nvar, 1L)
  varsel_star <- varsel
  varsel_star[change_ind] <- !varsel_star[change_ind]
  change_name <- names(regime_object)[change_ind]
  
  if( varsel[change_ind] ) {
    
    # Move to a lower dimension
    Theta_drop <- c(beta_new[change_name], 
                    regime_object[[change_name]]$DLF_parameters)
    T_drop <- Theta_drop
    T_drop <- log(-Theta_drop[3])
    beta_drop <- beta_new[setdiff(names(beta_new), change_name)]
    
    X_mat_star <- make_design_mat( regime_object, n_y, varsel_star )
    Xb_star <- X_mat_star %*% beta_drop
    ll_star <- -0.5*(crossprod(Xb_star) - z_pn%*%Xb_star)
    
    lDLF_prior_star <- 0
    for (k in regime_object[varsel_star]) {
      DLF_pars <- k$DLF_parameters
      lDLF_prior_star <- sum(lDLF_prior_star, dnorm(DLF_pars[1], sd = 10, log = TRUE), 
                             log(dtruncnorm(DLF_pars[2], b=0, sd = 10)))
    }
    
    lpost_star <- sum( ll_star,
                      dnorm(beta_drop, 0, 10, log = TRUE),
                      sum(dbinom(varsel_star, 1, varsel_prior, log = TRUE)),
                      lDLF_prior_star )
    
    lproposal <- sum(dnorm(T_drop, 0, 0.2, log = TRUE))
    accept_prob <- sum( lpost_star, -lpost, lproposal, -Theta_drop[3] )
    
    if (accept_prob >= log( runif(1) )) {
      varsel <<- vsres[i, ] <<- varsel_star
      beta_new <<- beta_drop
      X_mat <<- X_mat_star
      ll <<-ll_star
    }
    else {
      vsres[i, ] <<- vsres[i - 1L, ]
    }
  }
  else {
    # Move to a higher dimension
    beta_star_new <- rnorm(1, 0, 0.2)
    beta_star <- ins(beta_new, beta_star_new, change_ind)
    DLF_star <- rnorm(2, 0, 0.2)
    DLF_star_trans <- DLF_star
    DLF_star_trans[2] <- ljacob <- -exp(DLF_star[2])
    
    regime_object[[change_name]]$DLF_parameters <- DLF_star_trans
    X_mat_star <- make_design_mat( regime_object, n_y, varsel_star )
    Xb_star <- X_mat_star %*% beta_star 
    ll_star <- -0.5*(crossprod(Xb_star) - z_pn%*%Xb_star)
    
    lDLF_prior_star <- 0
    for (k in regime_object[varsel_star]) {
      DLF_pars <- k$DLF_parameters
      lDLF_prior_star <- sum(lDLF_prior_star, dnorm(DLF_pars[1], sd = 10, log = TRUE), 
                             log(dtruncnorm(DLF_pars[2], b=0, sd = 10)))
    }
    
    lpost_star <- sum( ll_star,
                      dnorm(beta_star, 0, 10, log = TRUE),
                      sum(dbinom(varsel_star, 1, varsel_prior, log = TRUE)),
                      lDLF_prior_star )
    
    lproposal <- sum(dnorm(c(beta_star_new, DLF_star), 0, 0.2, log = TRUE))
    accept_prob <- sum( lpost_star, -lpost,
                        -lproposal, ljacob )
    
    if (accept_prob >= log( runif(1) )) {
      varsel <<- vsres[i, ] <<- varsel_star
      beta_new <<- beta_star
      names(beta_new) <<- c( "Int", names(regime_object[varsel]) )
      ll <<- ll_star
      X_mat <<- X_mat_star
      regime_object[[change_ind]]$DLF_parameters <<- DLF_star_trans
    }
    else {
      vsres[i, ] <<- vsres[i - 1L, ]
    }
  }
  
}
