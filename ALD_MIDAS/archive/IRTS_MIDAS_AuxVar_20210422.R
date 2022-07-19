###### Auxiliary Variable Binary Inference
###### Daniel Dempsey

IRTS_MIDAS_AuxVar <- function(formula, data, response_dist = "ald", quantile = 0.5, prior, 
                              beta_start, DLF_start, incl_DLF_pars = TRUE, varsel = FALSE, 
                              MCMC_length = 10000L, ...) {
  
  #### For the moment, NEED to specify DLF prior and starting values to work
  
  ### Initialise environment to be passed into model.frame
  if (missing(data))
    Zenv <- environment()
  else
    Zenv <- as.environment(data)
  parent.env(Zenv) <- environment(formula)
  
  ### Check response_dist is valid
  valid_settings <- c("normal", "logistic", "ald", "error")
  chosen_setting <- pmatch(response_dist, valid_settings, 4)
  response_dist <- valid_settings[chosen_setting]
  if ( response_dist == "error" )
    stop("Supported values for response_dist are 'normal', 'logistic' and 'ald'.")
  
  ### Check if DLF parameter inference is desired
  if ( !incl_DLF_pars )
    response_dist <- paste0( response_dist, '_fixedDLF' )
  
  ### Check if variable selection is desired
  if ( varsel )
    response_dist <- paste0( response_dist, '_vs' ) 
  
  ### Unpack the formula
  form <- midas_formula_unpack(formula, Zenv)
  
  ### Check prior
  if ( missing(prior) ) 
    prior <- get(paste0("defaultPrior_", response_dist))(form, Zenv)
  
  ### Check starting values
  # Beta values
  n_var <- length( prior$beta0 )
  if ( missing(beta_start) )
    beta_start <- prior$beta0
  if ( length(start) < n_var ) {
    warning("Too few values in beta_start. Using mean of prior distribution instead.")
    beta_start <- prior$beta0
  }
  if ( length(start) > n_var ) {
    warning("Too many values in beta_start. Truncating to suitable length.")
    beta_start <- beta_start[1:n_var]
  }
  start <- list(beta = start_beta)
  
  # DLF parameter values
  # What values would be good????
  if ( incl_DLF_pars ) {
    if ( missing(DLF_start) )
      DLF_start <- 0
    start$DLF <- DLF_start
  }
  
  ### Run the MCMC routine
  routine_fun <- get(paste0("MIDAS_MCMC_", response_dist))
  routine_fun(formula = form, start = start, prior = prior, quantile = quantile, MCMC_length = MCMC_length, ...)
  
}

if (FALSE) {
  set.seed(123)
  binom_innov <- function(n, x) {
    fam <- binomial()
    rbinom(n, prob = fam$linkinv(x), size = 1)
  }
  n_vars <- 9
  sim <- irts_midas_sim(innov_fun = binom_innov, n_num = n_vars, beta = c(-2.5, 0, 1, -1, 0, 2, 0, -2, 0, 0.5))
  
  yy <- na.omit( sim$Data$y ) 
  qq <- 0.5
  
  mrtest <- midas_regimes(value = sim$Data[, seq(1, n_vars)], time = sim$Data$INDEX, DLF_parameters = list(c(0.5, -0.1)))
  rtest <- reponse_regime(value = sim$Data$y, time = sim$Data$INDEX)
  form <- rtest ~ mrtest
  
  ### Bayes QR fit
  Zenv <- as.environment(sim$Data)
  parent.env(Zenv) <- environment()
  form2 <- midas_formula_unpack(form, Zenv)
  library(bayesQR)
  tic()
  resqr <- bayesQR(formula = form2, data = sim$Data, quantile = qq, normal.approx = FALSE, ndraw = 10000L)
  toc()
  
  # Without variable selection
  tic()
  res1 <- IRTS_MIDAS_AuxVar(form, data = sim$Data, quantile = qq, response_dist = "ald", varsel = FALSE)
  toc()
  
  # With variable selection
  tic()
  #myprior <- list(beta0 = rep(0, 4), V0 = diag(100, 4), vars = c(1, 0.01, 0.99, 0.99))
  res2 <- IRTS_MIDAS_AuxVar(form, data = sim$Data, quantile = 0.5, response_dist = "ald", varsel = TRUE, MCMC_length = 10000)
  toc()
  
  res3 <- IRTS_MIDAS_AuxVar(form, data = sim$Data, quantile = qq, response_dist = "ald", varsel = TRUE, MCMC_length = 10000)
}

