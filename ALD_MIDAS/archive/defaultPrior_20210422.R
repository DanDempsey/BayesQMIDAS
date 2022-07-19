###### Default Priors for IRTS MIDAS Parameters
###### Daniel Dempsey

defaultPrior <- function(formula, data, incl_DLF_pars, varsel) {
  
  startm <- model.matrix( formula, data )
  n_beta <- NCOL( startm )
  prior_beta <- list( beta0 = rep(0, n_beta), 
                      V0 = 100 * diag(n_beta) )
  
  prior <- list(beta = prior_beta)
  
  if ( varsel )
    prior$vars <- c(1, rep(0.5, n_beta - 1))
  
  # What default would be good????
  if ( incl_DLF_pars )
    prior$DLF_pars <- 0
  
  prior
  
}

defaultPrior_ald <- function(formula, data) {
  
  startm <- model.matrix( formula, data )
  n_beta <- NCOL( startm ) 
  list( beta0 = rep(0, n_beta), 
        V0 = 100 * diag(n_beta) )
  
}

defaultPrior_ald_vs <- function(formula, data) {
  
  startm <- model.matrix( formula, data )
  n_beta <- NCOL( startm ) 
  list( beta0 = rep(0, n_beta), 
        V0 = 100 * diag(n_beta), 
        vars = c(1, rep(0.5, n_beta - 1)) )
  
}
