##### Compute posterior predictive distribution for ALD_vs
##### Daniel Dempsey

##### Function that builds the MIDAS design matrix
run_design_mat <- function(regime_object, offset_only, vars_ind) {
  
  if ( !any(vars_ind) ) 
    return( offset_only )
  
  cbind( offset_only, midas_design_matrices(regime_object[vars_ind] ) )
  
}
  
### Function that creates design matrices for test set
create_test_X <- function(reg_obj, var_sel, DLF, n_y) {
  
  int <- matrix(1, nrow = n_y, ncol = 1)
  N <- ncol(var_sel)
  res <- vector('list', N) 
  for ( i in 1:N ) {
    reg_obj_update <- Map( function(x, y, i) { x$DLF_parameters <- y[i, ] ; 
                           return(x) }, 
                           x = reg_obj, y = DLF, MoreArgs = list(i = i) )
    res[[i]] <- run_design_mat(reg_obj_update, int, var_sel[, i])
  }
  
  res
  
}

### Posterior Predictive Distribution for ALD Variable Selection 
ppd_ALD_vs <- function(MCMC_res, new_form, burn_in = 1000L, thin = 10L, ...) {
 
  ### Extract list of parameters from model fit
  MCMC_length <- nrow( MCMC_res$betares )
  MCMC_extract <- seq(burn_in + 1L, MCMC_length, thin)
  gam <- as.data.frame( t( MCMC_res$vsres[MCMC_extract, ] ) )
  gam_int <- rbind( TRUE, gam )
  bet_full <- as.data.frame( t( MCMC_res$betares[MCMC_extract, ] ) )
  bet <- Map( '[', x = bet_full, y = gam_int )
  DLF <- lapply(MCMC_res$DLFres, function(x, y) { x[y, ] }, y = MCMC_extract )
  model_sel_prob <- MCMC_res$model_sel_prob[MCMC_extract]
  prior <- MCMC_res$prior
  
  ### Design matrices for training data
  regime_object_train <- MCMC_res$regime_object
  yold <- MCMC_res$response_vec
  X_mat_train <- create_test_X(regime_object_train, gam, DLF, length(yold)) 
  Xb_train <- do.call('cbind', Map( '%*%', x = X_mat_train, y = bet ))
  rm( X_mat_train, regime_object_train )
  
  ### Design matrices for test data
  Zenv <- environment()
  parent.env(Zenv) <- environment(new_form)
  form <- midas_formula_unpack(new_form, Zenv)
  regime_object <- Zenv$regime_object
  ynew <- Zenv$response_vector
  X_mat_test <- create_test_X(regime_object, gam, DLF, length(ynew)) 
  Xb_test <- do.call('cbind', Map( '%*%', x = X_mat_test, y = bet )) 
  rm( X_mat_test, regime_object )
  
  ### Compute Bernoulli likelihood parameters
  train_p <- pALD(Xb_train, p = MCMC_res$quantile)
  test_p <- pALD(Xb_test, p = MCMC_res$quantile)
  
  ### Compute likelihood
  bern_train <- colSums( dbinom(yold, 1, train_p, log = TRUE) ) 
  
  ### Compute priors
  # beta given gamma
  prior_bet_list <- lapply(gam_int, function(x) { prior$beta0[x] })
  prior_sig_list <- lapply(gam_int, function(x) { as.matrix(prior$V0[x, x]) })
  prior_beta <- mapply( dmvnorm, x = bet, mean = prior_bet_list, sigma = prior_sig_list, 
                        MoreArgs = list(log = TRUE, checkSymmetry = FALSE ) )
  
  # theta given gamma
  prior_theta_list <- lapply( DLF, function(x) { apply(x, 1, all_DLF_prior, prior = prior) } )
  prior_theta_adjust <- Map( function(x, y) { as.numeric(x) * y }, 
                             x = as.data.frame(t(gam)), y = prior_theta_list )
  prior_theta <- rowSums( as.data.frame(prior_theta_adjust) )
  
  # gamma
  prior_gamma_list <- Map( dbinom, x = gam, prob = model_sel_prob, 
                           MoreArgs = list(size = 1, log = TRUE) )
  prior_gamma <- sapply(prior_gamma_list, sum)
  
  hyperprior_gamma <- sum( dbeta( model_sel_prob, prior$model_selection[1], 
                                  prior$model_selection[2] ) )
  
  ### Compute weights
  l_prop_w <- bern_train + prior_beta + prior_theta + prior_gamma + hyperprior_gamma
  w <- exp( l_prop_w - matrixStats::logSumExp(l_prop_w) )
  
  ### Compute posterior predictive probabilities and ROC
  ppp <- as.numeric( test_p %*% w )
  ROC <- pROC::roc( ynew, ppp, ... )
  
  ### Return results
  list( predicted = ppp, actual = ynew, ROC = ROC )
  
}
