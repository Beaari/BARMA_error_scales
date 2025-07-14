# =============================================================================
# Code for performing \betaARMA 
# =============================================================================
#
# Description:
#   Fits a Beta AutoRegressive Moving Average (BARMA) model for bounded time 
#   series data (in (0, 1)). Supports different link functions and computes 
#   forecasts, parameter estimates, Fisher information matrix, and model fit stats.
#
# Parameters:
#   y        : Time series object of proportions/rates in (0, 1)
#   ar       : Vector of AR lags (e.g., c(1, 2))
#   ma       : Vector of MA lags (e.g., c(1))
#   link     : Link function to use ("logit", "probit", "cloglog", or "loglog")
#   scale    : Error scale ("original" or "predictor")
#   h        : Forecast horizon (number of steps ahead)
#   X        : Optional covariate matrix (currently not used in core)
#   X_hat    : Optional covariate matrix for forecasting
#   ...      : Additional arguments (currently unused)
#
# Returns:
#   A list with model estimates, fitted values, forecasts, diagnostic matrices,
#   and information criteria.

barma <- function(y, 
                  ar = NA, ma = NA, 
                  link = "logit",
                  scale = "original",
                  h = 6, 
                  X = NA, X_hat = NA,...) {
  
  # Link Function, Inverse, and mu.eta  ----------------------------------------
  
  # Link structure
  link_structure <- link_choice(link)
  
  # Checking y structure                ----------------------------------------
  
  if (min(y) <= 0 || max(y) >= 1) stop("OUT OF RANGE (0,1)!")
  
  if (is.ts(y) == TRUE) {
    freq <- frequency(y)
  } else {
    stop("Data must be a time-series object")
  }
  
  
  # Parameters                          ----------------------------------------
  
  # Global parameters
  parameters <- list(
    link = link,
    scale = scale,
    
    ar = ar,
    ma = ma
  )
  
  # Adding additional calculable parameters to global parameter list
  parameters$linkfun <- link_structure$linkfun
  parameters$linkinv <- link_structure$linkinv
  parameters$mu.eta <- link_structure$mu.eta
  
  parameters$p = max(ar)
  parameters$q = max(ma)
  parameters$n = length(y)
  
  parameters$p_length = length(ar)
  parameters$q_length = length(ma)
  parameters$m = max(parameters$p, parameters$q, na.rm = T)
  
  parameters$y <-  y
  parameters$ynew <- link_structure$linkfun(y)
  parameters$ystar <- log(y / (1 - y))
  parameters$y_prev <- c(rep(NA, (parameters$n + h)))
  
  # Fetching parameters to be used in the rest of the code (excluding functions)
  n <- parameters$n
  m <- parameters$m
  
  p_length <- parameters$p_length
  q_length <- parameters$q_length
  
  ynew <- parameters$ynew
  y_prev <- parameters$y_prev
  
  linkfun <- parameters$linkfun
  linkinv <- parameters$linkinv
  mu.eta <- parameters$mu.eta
  
  # Renaming AR and MA parameters
  if (any(is.na(parameters$ar)) == FALSE) names_varphi <- paste("varphi", 
                                                                parameters$ar,
                                                                sep = "")
  
  if (any(is.na(parameters$ma)) == FALSE) names_theta  <- paste("theta", 
                                                                parameters$ma,
                                                                sep = "")
  
  
  # Start Values                        ----------------------------------------
  
  result <- list()
  
  result$n_obs <- n
  
  start_values <- startValues(y,
                              link = link,
                              parameters = parameters,
                              ar = ar,
                              ma = ma,
                              X = X)
  
  
  # Optimization                        ----------------------------------------
  
  names_par_arma <- c("alpha", names_varphi, names_theta, "phi")
  
  opt_arma <- stats::optim(
    par     = start_values,
    fn      = loglik_arma,
    gr      = score_arma, 
    method  = "BFGS",
    
    control = list(
      maxit = 1000,
      reltol = 1e-12
    ),
    parameters = parameters
  )
  
  if (opt_arma$conv != 0) {
    warning("FUNCTION DID NOT CONVERGE!")
  }
  
  # Convergence Status
  result$conv   <- opt_arma$convergence
  result$counts <- as.numeric(opt_arma$counts[1])
  result$opt    <- opt_arma
  
  # Estimates
  coef_arma        <- opt_arma$par[1:(p_length + q_length + 2)]
  names(coef_arma) <- c("alpha", names_varphi, names_theta, "phi")
  
  # Output Estimates
  result$coef <- coef_arma
  
  # Log Likelihood Value 
  result$loglik <- (-1) * opt_arma$value
  
  
  # Information Fisher Matrix           ----------------------------------------
  
  alpha  <- coef_arma[1]
  varphi <- coef_arma[2:(p_length + 1)]
  theta  <- coef_arma[(p_length + 2):(p_length + q_length + 1)]
  phi    <- coef_arma[p_length + q_length + 2]
  
  result$alpha  <- alpha
  result$varphi <- varphi
  result$theta  <- theta
  result$phi    <- phi
  
  # Fisher Information Matrix, ARMA
  output_inf_matrix_arma <- inf_matrix_arma(y = y,
                                            ar = ar,
                                            ma = ma,
                                            
                                            parameters = parameters,
                                            
                                            alpha = result$alpha, 
                                            varphi = result$varphi, 
                                            theta = result$theta,
                                            phi = result$phi,
                                            
                                            link = link)
  
  fisher_info_mat <- output_inf_matrix_arma$fisher_info_mat
  
  # output
  result$fisher_info_mat <- fisher_info_mat
  result$muhat <- output_inf_matrix_arma$muhat
  result$fitted <- output_inf_matrix_arma$fitted
  result$etahat <- output_inf_matrix_arma$etahat
  
  
  # Forecasting                         ----------------------------------------
  
  errorhat    <- output_inf_matrix_arma$errorhat
  ynew_prev <- c(ynew, rep(NA, h))
  y_prev[1:n] <- result$fitted
  
  for (i in 1:h) {
    
    ynew_prev[n + i] <- alpha + 
      (varphi %*% ynew_prev[n + i - ar]) + 
      theta %*% errorhat[n + i - ma]
    
    y_prev[n + i] <- linkinv(ynew_prev[n + i])
    
    errorhat[n + i] <- 0 
  }
  
  # Updating y_prev in parameters
  parameters$y_prev <- y_prev
  
  # Final Values                        ----------------------------------------
  
  output <- final_output(parameters, result, X, start_values, fisher_info_mat, h)
  
  return(output)
}