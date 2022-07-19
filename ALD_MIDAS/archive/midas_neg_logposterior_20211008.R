##### Computes the negative log posterior at given co-ordinates
##### Daniel Dempsey

midas_neg_logposterior <- function(pars, pinds, formula, family, prior, lp) {
  
  ### Extract formula environment
  Zenv <- environment(formula)
  
  ### Set the parameters parameters
  par_list <- split(pars, pinds)
  Zenv$regime_object <- mapply("[[<-", x = Zenv$regime_object, i = "DLF_parameters", value = par_list[-1], SIMPLIFY = FALSE)
  
  ### Compute the linear predictor
  WX <- model.matrix( formula, Zenv )
  WXb <- WX %*% par_list[[1]]
  eta <- family$linkinv(WXb)
  
  ### Compute the prior
  log_prior <- prior(pars, pinds)
  if (!lp)
    log_prior <- log(log_prior)
  
  ### Return log posterior density
  sum( family$dev.resids(Zenv$yv, eta, 1L) ) - 2*log_prior
  
}

if (FALSE) {
  
  set.seed(123)
  test <- irts_midas_sim(n_num = 3, n_fac = 0, n_levs = 3, y_n = 200, x_n = 600, x_lev_effect = lev_effect,
                         x_num_effect = 4, theta = c(0.5, -0.1), tau = 10, family = "gaussian")
  X <- test[-NCOL(test)]
  y <- na.omit(test[c(1, NCOL(test))])
  yt <- y$time
  yt2 <- yt
  y <- y$y
  Xt <- irts_time_matrix(X)
  mat_list <- irts_var_weights(Xt, yt, tau = 10, time_weights = 1)
  X <- X[-1]
  SSDM <- irts_design_matrix(X)
  ##
  #Loss(theta, mat_list, X, SSDM, y, indvec, theta_equal = FALSE, DLF = irts_nealmon)
  
}

