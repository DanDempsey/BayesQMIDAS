##### Compute posterior predictive distribution for ALD_vs
##### Daniel Dempsey
  
### Posterior Predictive Distribution for ALD Variable Selection 
ppd_ALD_vs <- function(MCMC_res, xnew, ynew = NULL, 
                       burn_in = 1000L, thin = 10L, ...) {
  
  MCMC_length <- nrow( MCMC_res$betares )
  MCMC_extract <- seq(burn_in + 1L, MCMC_length, thin)
  bet <- t( MCMC_res$betares[MCMC_extract, ] )
  gam <- as.data.frame( t( MCMC_res$vsres[MCMC_extract, ] ) )
  
  # Extract training and test set
  XB_train <- MCMC_res$design_mat %*% bet
  XB_test <- as.matrix(cbind(1, xnew)) %*% bet
  
  # Compute Bernoulli likelihood parameters
  train_p <- pALD(XB_train, p = MCMC_res$quantile)
  test_p <- pALD(XB_test, p = MCMC_res$quantile)
  
  # Compute likelihood
  bern_train <- colSums( dbinom(MCMC_res$response_vec, 1, train_p, log = TRUE) ) 
  
  # Compute beta prior given gamma
  prior_bet_list <- lapply(gam, function(x) { MCMC_res$prior$beta0[x] })
  prior_sig_list <- lapply(gam, function(x) { MCMC_res$prior$V0[x, x] })
  bet_draw_list <- Map(function(x, y) { x[y] }, x = as.data.frame(bet), y = gam)
  
  prior_vals <- mapply( dmvnorm, x = bet_draw_list, 
                        mean = prior_bet_list, 
                        sigma = prior_sig_list,
                        MoreArgs = list(log = TRUE) )
  
  # Compute weights
  l_prop_w <- bern_train + prior_vals + sum(log(MCMC_res$prior$vars))
  w <- exp( l_prop_w - matrixStats::logSumExp(l_prop_w) )
  
  # Compute posterior predictive probabilities
  ppp <- as.numeric( test_p %*% w )
  if (is.null(ynew)) { return(ppp) }
  ROC <- pROC::roc(ynew, ppp, ...)
  
  list(predicted = ppp, actual = ynew, ROC = ROC)
  
}
