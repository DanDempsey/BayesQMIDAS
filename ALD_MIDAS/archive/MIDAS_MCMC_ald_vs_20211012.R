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
make_design_mat <- function(regime_object, offset_only, vars_ind) {
  
  if ( !any(vars_ind) ) 
    return( offset_only )
  
  cbind( offset_only, midas_design_matrices(regime_object[vars_ind] ) )
  
}

##### Function that computes DLF prior
all_DLF_prior <- function(x, prior) {
  DLF_pars <- prior$DLF_pars
  sum( dnorm(x[1], mean = DLF_pars[[1]][1], sd = DLF_pars[[1]][2], log = TRUE), 
       dgamma(-x[2], shape = DLF_pars[[2]][1], rate = DLF_pars[[2]][2], log = TRUE) )
}

##### Function that computes covariance matrix
compute_covar <- function(x, y, z, i, j) {
  x[, 2] <- log(-x[, 2])
  z * cov(x[y[1:(i-1), j], ])
}

##### Main function that performs inference via MCMC
MIDAS_MCMC_ald_vs <- function(formula, start, prior, quantile, MCMC_length) {
  
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
  X_mat <- make_design_mat( regime_object, int, varsel )
  nu <- rep( 1, n_y )
  Ni <- 1 / ( omega * nu )
  z <- rep(0, n_y)
  z_pn <- z - (psi * nu) 
  
  # Beta initialization
  betares <- matrix( 0, ncol = nvar+1, nrow = MCMC_length )
  vsres <- matrix( FALSE, ncol = nvar, nrow = MCMC_length )
  model_sel_prob <- matrix( 0.5, ncol = nvar, nrow = MCMC_length )
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
  birth_oppurtunity <- birth <- death_oppurtunity <- death <- accept <- rep( 0, nvar )
  sigmaDLF_start <- rep( 0.1, nvar )
  covar_sd <- (2.38^2)/2
  DLF_propose_covar <- replicate( nvar, covar_sd * matrix(c(6, 0.5, 0.5, 0.05), ncol = 2), simplify = FALSE )
  names(accept) <-  names(birth_oppurtunity) <- 
    names(birth) <- names(death_oppurtunity) <- names(death) <- 
    names(sigmaDLF_start) <- varnames 
  rm(var_extract, DLFres_cmd, DLF_1, varnames)
  
  ### Main Loop
  for ( i in 2:MCMC_length ) {
    
    if ( i == 100 ) { browser() }
    ### Update covariate indicator
    V0i <- V0i_full[varsel_int, varsel_int]
    V0ib0 <- V0ib0_full[varsel_int]
    XtNi <- t( X_mat * Ni )
    
    V_posti <- V0i + XtNi%*%X_mat
    V_post <- chol2inv( chol(V_posti) )
    B_post <- V_post%*%( V0ib0 + XtNi%*%z_pn )
    
    if ( i%%10 != 0 ) { 
      
      vsres[i, ] <- vsres[i-1, ]
      model_sel_prob[i, ] <- model_sel_prob[i-1, ]
      
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
        
        sigmaDLF_j <- sigmaDLF_start[change_name]
        
        # Move to a higher dimension
        Theta_star <- rnorm(2, 0, sigmaDLF_j)
        lproposal <- dnorm(Theta_star, 0, sigmaDLF_j, log = TRUE)
        ljacob <- Theta_star[2]
        Theta_star[2] <- -exp(Theta_star[2])
        regime_object[[change_name]]$DLF_parameters <- Theta_star
        lDLFprior <- all_DLF_prior( Theta_star, prior )
        DLF_component <- sum( -lproposal, lDLFprior, ljacob )
        birth_oppurtunity[change_name] <- birth_oppurtunity[change_name] + 1
        
      }
      
      X_mat_star <- make_design_mat( regime_object, int, varsel_star )
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
      varsel_lprior <- dbinom(varsel, 1, model_sel_prob[i-1, ], log = TRUE)
      varsel_lprior_star <- dbinom(varsel_star, 1, model_sel_prob[i-1, ], log = TRUE)
      
      lkernel <- crossprod(B_post, V_posti)%*%B_post/2
      lkernel_star <- crossprod(B_post_star, V_posti_star)%*%B_post_star/2
      ldenom <- sum( ldet_V_post, ldet_V0i_star, lkernel, varsel_lprior )
      lnum <- sum( ldet_V_post_star, ldet_V0i, lkernel_star, 
                   varsel_lprior_star, DLF_component ) 
      
      if ( (lnum - ldenom) > log(runif(1)) ) {
        if( !varsel[change_ind] ) {
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
      
      ### Update the prior success probability p
      model_sel_prob[i, ] <- rbeta(nvar, 1 + varsel, 2 - varsel)
        
    }
    
    ### Update ALD parameters
    # Update beta
    bet <- betares[i, varsel_int] <- t( rmvnorm( 1, mean = B_post, sigma = V_post, checkSymmetry = FALSE ) )
    
    # Update z
    Xb <- X_mat %*% bet
    z <- rTALD( n = n_y, rtrunc = rtrunc, mu = Xb, sigma = 1, p = quantile )
    
    # Update nu
    chi <- (z - Xb)^2 / omega
    nu <- 1/rinvgauss( n_y, mean = sqrt(chi/delta), shape = chi )
    
    ### Update DLF parameters
    z_pn <- z - (psi * nu)
    Ni <- 1 / ( omega * nu )
    z_pn_Ni <- z_pn * Ni
    ll <- z_pn_Ni%*%Xb - crossprod(Xb*Ni, Xb)/2
    
    # Loop over DLF updates
    for (j in which(varsel)) {
      
      # Acceptance Ratio Denominator
      DLF_last <- DLF_trans_last <- regime_object[[j]]$DLF_parameters
      l_r_denom <- ll + all_DLF_prior( DLF_last, prior )
      
      # Proposal
      if ( sum( vsres[1:(i-1), j] ) >= 100 ) { 
        DLF_propose_covar[[j]] <- compute_covar( DLFres[[j]][1:(i-1), ], vsres, covar_sd, i, j ) 
      }
      DLF_trans_last[2] <- log( -DLF_trans_last[2] )
      DLF_star <- t( rmvnorm( 1, mean = DLF_trans_last, sigma = DLF_propose_covar[[j]] ) )
      ljacob <- DLF_star[2] - DLF_trans_last[2]
      DLF_star[2] <- -exp( DLF_star[2] )
      
      # Acceptance Ratio Numerator
      regime_object[[j]]$DLF_parameters <- DLF_star
      X_mat_star <- cbind( int, midas_design_matrices(regime_object[varsel]) )
      Xb_star <- X_mat_star %*% bet
      ll_star <- z_pn_Ni%*%Xb_star - crossprod(Xb_star*Ni, Xb_star)/2
      l_r_num <- ll_star + all_DLF_prior( DLF_star, prior ) + ljacob
      
      if ( (l_r_num - l_r_denom) > log(runif(1)) ) {
        DLFres[[j]][i, ] <- DLF_star
        accept[j] <- accept[j] + 1
        ll <- ll_star
        X_mat <- X_mat_star
      } else {
        DLFres[[j]][i, ] <- regime_object[[j]]$DLF_parameters <- DLF_last
      }
    }
    
    # Print progress
    if ( !i%%500 )
      cat( paste0("Current iteration: ", i, "\n") )

  }
  
  list(betares = betares, DLFres = DLFres, vsres = vsres,
       model_sel_prob = model_sel_prob,
       quantile = quantile, response_vec = Zenv$response_vector, 
       prior = prior, accept = accept,
       sigmaDLF_start = sigmaDLF_start, 
       birth_oppurtunity = birth_oppurtunity,
       death_oppurtunity = death_oppurtunity,
       birth = birth, death = death, 
       regime_object = regime_object)
  
}
