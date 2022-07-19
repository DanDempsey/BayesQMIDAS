##### Infers the Maximum A Posteriori of the IRTS-MIDAS
##### Daniel Dempsey

irts_midas_MAP <- function(formula, data, group = NULL, start = NULL, family = "gaussian", prior = nealmon_prior, 
                           logprior = TRUE, n_cores = 1L, Ofunction = "optim", ...) {
  
  ### Initialise environment to be passed into model.frame
  if (missing(data))
    Zenv <- new.env(parent = environment())
  else {
    Zenv <- as.environment(data)
    parent.env(Zenv) <- environment()
  }
  
  ### Check that the family parameter is valid
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    stop("Supplied family is not recognized.")
  }
  
  ### Check that Ofunction parameter is valid
  if ( is.function(Ofunction) )
    Ofunction <- as.character( substitute(Ofunction) )
  if( !is.character(Ofunction) )
    stop("Ofunction must be either a function or string.")
  if( !(Ofunction %in% c("optim", "spg", "optimx", "dry_run")) ) 
    stop("Supplied Ofunction is not in the supported optimization functions list.")
  
  ### Unpack the formula
  form <- midas_formula_unpack(formula, Zenv)
  
  ### Assign starting values and prepare parameters
  MIDAS_start(form, Zenv)
  
  ### Gather all relevant info and run the fitting function
  Ofunction_pars <- list(Ofunction = Ofunction, par = Zenv$pars, fn = midas_neg_logposterior, family = family, 
                         prior = prior, lp = logprior, pinds = Zenv$pinds, formula = form, ...)
  invisible( list2env(list(Ofunction = Ofunction, n_cores = n_cores), envir = Zenv) )
  irts_midas_MAP.fit(Ofunction_pars)
  
}

if(FALSE){
  
  binom_innov <- function(n, x) {
    fam <- binomial()
    rbinom(n, prob = fam$linkinv(x), size = 1)
  }
  
  set.seed(123)
  test <- irts_midas_sim(n_y = 1500, n_num = 1, n_fac = 1, beta = c(0, 2, -2), innov_fun = binom_innov)$Data
  
  mrtest <- midas_regimes(value = test[1:2], time = test$INDEX, DLF_parameters = list(c(0, -0.1)))
  form <- y ~ mrtest | INDEX
  
  #tic()
  fit <- irts_midas_MAP(form, data = test, n_cores = 1, family = "binomial", method = "BFGS", hessian = TRUE)
  #toc()
  
  ## Without hessian:
  # 1 cores: 8s
  # 2 cores: 15s
  # 3 cores: 15s
  # 4 cores: 14s
  
  ## With hessian:
  # 1 cores: 13s
  # 2 cores: 24s
  # 3 cores: 24s
  # 4 cores: 24s
  
  ### Non-DLF term?
  ndlf <- rnorm(test$y %>% na.omit() %>% length())
  ndlf2 <- rnorm(test$y %>% na.omit() %>% length())
  ndlf3 <- rnorm(test$y %>% na.omit() %>% length())
  
  form <- y ~ ndlf + ndlf2 + ndlf3 + mrtest | INDEX
  fit <- irts_midas_MAP(form, data = test, n_cores = 1, family = "binomial", method = "BFGS", hessian = TRUE)
  
}


