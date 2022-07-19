####### ALD MIDAS Regression with Reversible-Jump Variable Selection 
####### Daniel Dempsey

##### Function that inserts a value into a specific location in a vector
ins <- function(a, to.insert, pos) {

  pos <- pos + 1
  if (pos <= 1)
    return( c(to.insert, a[seq(pos[1], length(a))]) )
  if (pos >= length(a)+1)
    return( c(a, to.insert)  )
  
  c(a[seq(pos[1]-1)], to.insert, a[seq(pos[1], length(a))])
  
}

##### Function that builds the MIDAS design matrix
make_design_mat <- function(regime_object, offset_only, vars_ind, n_cores) {
  
  if ( !any(vars_ind) ) 
    return( offset_only )
  
  cbind( offset_only, midas_design_matrices(regime_object[vars_ind], n_cores) )
  
}

##### Function that computes DLF prior
all_DLF_prior <- function(x, prior) {
  DLF_pars <- prior$DLF_pars
  sum( dnorm(x[1], mean = DLF_pars[[1]][1], sd = DLF_pars[[1]][2], log = TRUE), 
       dgamma(-x[2], shape = DLF_pars[[2]][1], rate = DLF_pars[[2]][2], log = TRUE) )
}

##### Function that extracts DLF parameters 1 and 2
DLF_par_extract <- function( regime_object, varsel ) {
  reg_select <- regime_object[varsel]
  if ( length(reg_select) == 0 ) {
    return( NA )
  }
  sapply( reg_select, function(x) { x$DLF_parameters } )
}

##### Function that computes separate DLF priors
DLF_prior <- function( x, type, prior ) {
  if (type == 1) {
    DLF_pars <- prior$DLF_pars[[1]]
    return( sum( dnorm(x, mean = DLF_pars[1], sd = DLF_pars[2], log = TRUE) ) )
  }
  DLF_pars <- prior$DLF_pars[[2]]
  sum( dgamma(-x, shape = DLF_pars[1], rate = DLF_pars[2]) )
}

##### Function that sets DLF values to regime object
DLF_set <- function( pars, regime_object, type ) {
  regime_object$DLF_parameters[type] <- pars
  regime_object
}

##### Function that sets DLF results
DLFres_fun <- function(dat, row_ind, col_ind, x) { 
  dat[row_ind, col_ind] <- x 
  dat 
}

##### Main function that performs inference via MCMC
MIDAS_MCMC_ald_vs <- function(formula, start, prior, quantile, adapt, n_cores, 
                              MCMC_length, nburn) {
  
  ### Prepare data environment
  Zenv <- environment( formula )
  regime_object <- Zenv$regime_object
  nvar <- length(regime_object)
  
  ### MCMC setup (first iteration)
  rtrunc <- ifelse( Zenv$response_vector, TRUE, FALSE )
  n_y <- length(rtrunc)
  
  psi <- ( 1 - 2*quantile ) / ( quantile * (1 - quantile) )
  omega <- 2 / ( quantile * (1 - quantile) )
  delta <- 2 + ( ( psi^2 ) / omega )
  V0i_full <- chol2inv( chol(prior$V0) )
  V0i <- V0i_full[1, 1]
  V0ib0_full <- V0i_full %*% prior$beta0
  V0ib0 <- V0ib0_full[1, 1]
  
  int <- matrix(1, nrow = n_y, ncol = 1)
  varsel <- rep( FALSE, nvar )
  X_mat <- make_design_mat( regime_object, int, varsel, n_cores )
  nu <- rep( 1, n_y )
  Ni <- 1 / ( omega * nu )
  z <- rep(0, n_y)
  z_pn <- z - (psi * nu) 
  
  # Beta initialization
  betares <- matrix(0, ncol = nvar+1, nrow = MCMC_length)
  vsres <- matrix(FALSE, ncol = nvar, nrow = MCMC_length)
  varnames <- names(Zenv$regime_object)
  colnames(betares) <- c( "Int", varnames )
  colnames(vsres) <- varnames
  betares[1, ] <- start
  varsel <- vsres[1, ]
  varsel_int <- c( TRUE, varsel )
  
  # DLF parameter initialization
  var_extract <- paste0('regime_object$', varnames, '$DLF_parameters')
  DLFres_cmd <- Map( parse, text = var_extract )
  DLF_1 <- lapply(DLFres_cmd, eval, envir = Zenv)
  names(DLF_1) <- varnames
  DLFres <- lapply( DLF_1, function(x) { matrix(rep(x, each = MCMC_length), nrow = MCMC_length) } )
  
  var_inds <- 1:nvar
  accept <- runs <- rep( 1, 2 )
  accept_start <- runs_start <- rep( 1, nvar )
  accept_after_burn_in <- runs_after_burn_in <- rep( 0, 2 )
  birth_oppurtunity <- birth <- death_oppurtunity <- death <- rep( 0, nvar )
  sigmaDLF <- rep( 0.1, 2 )
  sigmaDLF_start <- rep( 0.1, nvar )
  names(sigmaDLF) <- names(runs) <- names(accept) <- 
    names(accept_after_burn_in) <- names(runs_after_burn_in) <- c('DLF1', 'DLF2')
  names(birth_oppurtunity) <- names(birth) <- names(death_oppurtunity) <- 
    names(death) <- names(runs_start) <- names(sigmaDLF_start) <- 
    names(accept_start) <- varnames 
  rm(var_extract, DLFres_cmd, DLF_1, varnames)
  
  ### Main Loop
  for ( i in 2:MCMC_length ) {
    
    #if ( i == 510 ) { browser() }
    ### Update covariate indicator
    V0i <- V0i_full[varsel_int, varsel_int]
    V0ib0 <- V0ib0_full[varsel_int]
    XtNi <- t( X_mat * Ni )
    
    V_posti <- V0i + XtNi%*%X_mat
    V_post <- chol2inv( chol(V_posti) )
    B_post <- V_post%*%( V0ib0 + XtNi%*%z_pn )
    
    if ( i%%10 != 0 ) { 
      
      vsres[i, ] <- vsres[i-1, ]
      
    } else {
      
      ### Propose dimension change
      change_ind <- sample( var_inds, 1 )
      varsel_star <- varsel
      varsel_star[change_ind] <- !varsel_star[change_ind]
      varsel_star_int <- c( TRUE, varsel_star )
      change_name <- names(regime_object)[change_ind]
      
      if( varsel[change_ind] ) {
        
        # Move to a lower dimension
        sigmaDLF_j <- sigmaDLF_start[change_name]
        Theta_drop <- regime_object[[change_name]]$DLF_parameters
        lDLFprior <- all_DLF_prior( Theta_drop, prior )
        Theta_drop[2] <- log( -Theta_drop[2] )
        ljacob <- Theta_drop[2]
        lproposal <- dnorm(Theta_drop, sd = sigmaDLF_j, log = TRUE)
        DLF_component <- sum( lproposal, -lDLFprior, -ljacob )
        death_oppurtunity[change_name] <- death_oppurtunity[change_name] + 1
        
      } else {
        
        # Proposal variance adaptation
        #runs_j <- runs_start[change_name]
        #if( adapt & (i < nburn) & (runs_j%%100L == 0L) ){
        #  Delta <- 1/sqrt(runs_j)
        #  adapt_indicator <- ifelse( accept_start[change_name]/runs_j > .1, TRUE, FALSE ) # Don't expect to jump between models as often
        #  if (adapt_indicator)
        #    sigmaDLF_start[change_name] <- sigmaDLF_start[change_name] * exp( Delta )
        #  else
        #    sigmaDLF_start[change_name] <- sigmaDLF_start[change_name] * exp( -Delta )
        #}
        
        sigmaDLF_j <- sigmaDLF_start[change_name]
        #runs_start[change_name] <- runs_start[change_name] + 1L
        
        # Move to a higher dimension
        Theta_star <- rnorm(2, 0, sigmaDLF_j)
        lproposal <- dnorm(Theta_star, 0, sigmaDLF_j, log = TRUE)
        ljacob <- Theta_star[2]
        Theta_star[2] <- -exp(Theta_star[2])
        regime_object[[change_name]]$DLF_parameters <- Theta_star
        lDLFprior <- all_DLF_prior( Theta_star, prior )
        DLF_component <- sum( -lproposal, lDLFprior, ljacob )
        birth_oppurtunity[change_name] <- birth_oppurtunity[change_name] + 1L
        
      }
     
      X_mat_star <- make_design_mat( regime_object, int, varsel_star, n_cores )
      V0i_star <- V0i_full[varsel_star_int, varsel_star_int]
      V0ib0_star <- V0ib0_full[varsel_star_int]
      XtNi_star <- t( X_mat_star * Ni )
      
      V_posti_star <- V0i_star + XtNi_star%*%X_mat_star
      V_post_star <- chol2inv( chol(V_posti_star) ) 
      B_post_star <- V_post_star%*%( V0ib0_star + XtNi_star%*%z_pn )
      
      ldet_V_post <- sum( log(diag(chol(V_post))) ) 
      ldet_V_post_star <- sum( log(diag(chol(V_post_star))) )
      ldet_V0i <- sum( log(diag(chol(V0i))) )
      ldet_V0i_star <- sum( log(diag(chol(V0i_star))) )
      varsel_lprior <- dbinom(varsel_int, 1, prior$vars, log = TRUE)
      varsel_lprior_star <- dbinom(varsel_star_int, 1, prior$vars, log = TRUE)
      
      lkernel <- crossprod(B_post, V_posti)%*%B_post/2
      lkernel_star <- crossprod(B_post_star, V_posti_star)%*%B_post_star/2
      ldenom <- sum( ldet_V_post, ldet_V0i_star, lkernel, varsel_lprior )
      lnum <- sum( ldet_V_post_star, ldet_V0i, lkernel_star, 
                   varsel_lprior_star, DLF_component ) 
      
      if ( (lnum - ldenom) > log(runif(1L)) ) {
        if( !varsel[change_ind] ) {
          #accept_start[change_name] <- accept_start[change_name] + 1L
          birth[change_name] <- birth[change_name] + 1
        } else {
          DLFres[[change_name]][i, ] <- DLFres[[change_name]][1, ]
          death[change_name] <- death[change_name] + 1
        }
        vsres[i, ] <- varsel <- varsel_star
        varsel_int <- varsel_star_int
        X_mat <- X_mat_star
        V_post <- V_post_star
        B_post <- B_post_star
      } else {
        vsres[i, ] <- vsres[i-1, ]
      }
        
    }
    
    ### Update ALD parameters
    # Update beta
    bet <- betares[i, varsel_int] <- mvrnorm( mu = B_post, Sigma = V_post )
    
    # Update z
    Xb <- X_mat %*% bet
    z <- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
    
    # Update nu
    chi <- (z - Xb)^2 / omega
    for (k in 1:n_y) {
      nu[k] <- rgig( n = 1, lambda = 0.5, chi = chi[k], psi = delta )
    }
    
    ### Update DLF parameters
    z_pn <- z - (psi * nu)
    Ni <- 1 / ( omega * nu )
    z_pn_Ni <- z_pn * Ni
    ll <- z_pn_Ni%*%Xb - crossprod(Xb*Ni, Xb)/2
    
    # Loop over DLF updates
    par_mat <- DLF_par_extract( regime_object, varsel )
    if ( !all( is.na(par_mat) ) ) {
      for ( j in 1:2 ) {
        # Proposal variance adaptation
        if( adapt & (i < nburn) & (runs[j]%%200 == 0L) ){
          Delta <- 1/sqrt(runs[j])
          adapt_indicator <- ifelse( accept[j]/runs[j] > .234, TRUE, FALSE ) 
          if (adapt_indicator) {
            sigmaDLF[j] <- sigmaDLF[j] * exp( Delta )
          } else {
            sigmaDLF[j] <- sigmaDLF[j] * exp( -Delta )
          }
        }
        runs[j] <- runs[j] + 1
        if ( i >= nburn ) { runs_after_burn_in[j] <- runs_after_burn_in[j] + 1 }
        
        # Acceptance Ratio Denominator
        DLF_par <- par_mat[j, ]
        l_r_denom <- ll + DLF_prior( DLF_par, j, prior )
        
        # Proposal
        if ( j == 1 ) { DLF_par_trans <- DLF_par }
        else { DLF_par_trans <- log( -DLF_par ) }
        DLF_star <- rnorm( length(DLF_par), mean = DLF_par_trans, sd = sigmaDLF[j] )
        if ( j == 1 ) { ljacob <- 0 }
        else {
          ljacob <- sum( DLF_star - DLF_par_trans )
          DLF_star <- -exp( DLF_star )
        }
        
        # Acceptance Ratio Numerator
        regime_object[varsel] <- Map( DLF_set, pars = DLF_star, 
                                      regime_object[varsel], 
                                      MoreArgs = list(type = j) )
        X_mat_star <- cbind( int, midas_design_matrices(regime_object[varsel], n_cores) )
        Xb_star <- X_mat_star %*% bet
        ll_star <- z_pn_Ni%*%Xb_star - crossprod(Xb_star*Ni, Xb_star)/2
        l_r_num <- ll_star + DLF_prior( DLF_star, j, prior ) + ljacob
        
        # Accept or reject proposal
        if ( (l_r_num - l_r_denom) > log(runif(1)) ) {
          DLFres[varsel] <- Map( DLFres_fun, dat = DLFres[varsel], 
                                 row_ind = i, col_ind = j, x = DLF_star )
          accept[j] <- accept[j] + 1
          if ( i >= nburn ) { accept_after_burn_in[j] <- accept_after_burn_in[j] + 1 }
          ll <- ll_star
          X_mat <- X_mat_star
        } else {
          DLFres[varsel] <- Map( DLFres_fun, dat = DLFres[varsel],
                                 row_ind = i, col_ind = j, x = DLF_par )
          regime_object[varsel] <- Map( DLF_set, pars = DLF_par, 
                                      regime_object[varsel], 
                                      MoreArgs = list(type = j) )
        }
      }
    }
    
    # Print progress
    if ( !i%%500 )
      cat( paste0("Current iteration: ", i, "\n") )

  }
  
  list(betares = betares, DLFres = DLFres, vsres = vsres,
       quantile = quantile,
       response_vec = Zenv$response_vector, prior = prior,
       accept = accept, runs = runs, sigmaDLF = sigmaDLF,
       accept_after_burn_in = accept_after_burn_in, 
       runs_after_burn_in = runs_after_burn_in, 
       sigmaDLF_start = sigmaDLF_start, 
       birth_oppurtunity = birth_oppurtunity,
       death_oppurtunity = death_oppurtunity,
       birth = birth, death = death, 
       regime_object = regime_object)
  
}

