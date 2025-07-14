# Link Choice function               -------------------------------------------

# Description:
#   Selects and returns components of a link function for in BARMA.
#
#   This function provides the link function (`g(mu)`), its inverse (`g⁻¹(eta)`),
#   and the derivative of the inverse (`(g⁻¹)'(eta)`), depending on the specified
#   link type. It supports standard GLM link functions like "logit", "probit",
#   "cloglog", and "loglog".
#
# Parameters: 
#   link A character string indicating the link function to use.
#   Must be one of: "logit", "probit", "cloglog", "loglog".
#
# Output:
#  - link: The name of the selected link function.
#  - linkfun: The link function: g(mu).
#  - linkinv: The inverse link function: g⁻¹(eta).
#  - mu.eta: The derivative of the inverse link function: (g⁻¹)'(eta).

link_choice <- function(link){
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link") {
      linktemp <- eval(link)
    }
  }
  
  if (linktemp == "logit") {
    
    stats <- make.link("logit")
    
  } else if (linktemp == "probit") {
    
    stats <- make.link("probit")
    
  } else if (linktemp == "cloglog") {
    
    stats <- make.link("cloglog")
    
  } else if (linktemp == "loglog") {
    
    stats <- list()
    
    stats$linkfun <- function(mu) -log(-log(mu))
    stats$linkinv <- function(eta) exp(-exp(-eta))
    stats$mu.eta <- function(eta) exp(-exp(-eta)) * exp(-eta)
    
  } else {
    stop(
      paste0(
        "The '", linktemp, "' link function is not available.", "\n",
        "Available options are: ", "logit, loglog, cloglog, and probit."
      )
    )
  }
  
  link_structure <- structure(list(
    link = linktemp,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    mu.eta = stats$mu.eta
  ))
  
  return(link_structure)
}


# Derivatives per link               -------------------------------------------

# Description:
#   Computes the inverse of the derivative of the link function (1 / g'(mu))
#   for common link functions used in generalized linear models.
#
# Parameters:
#   mu   : Numeric vector of means (must be in (0, 1))
#   link : Character string specifying the link function. 
#          Supported options: "logit", "probit", "cloglog", "loglog"
#
# Output:
#   - A numeric value or vector with the computed 1 / g'(mu) for the specified link.
#   - Throws an error if the link function is not supported.

inv_g_prime <- function(mu, link) {
  if (link == "logit") {
    
    return(mu * (1 - mu))
    
  } else if (link == "probit") {
    
    z <- qnorm(mu)
    
    return(dnorm(z))
    
  } else if (link == "cloglog") {
    
    return((1 - mu) * log(1 - mu))
    
  } else if (link == "loglog") {
    
    return(mu * log(mu))
    
  } else {
    stop("Link function not supported")
  }
}


# Initial Values function            -------------------------------------------

# Description:
#   Computes starting values for the BARMA (Beta Autoregressive Moving Average) model 
#   parameters using a linear regression on the transformed response.
#   The function estimates initial values for alpha, AR coefficients (varphi),
#   MA coefficients (theta, initialized as zeros), and the precision parameter phi.
#
# Parameters:
#   y         : Numeric time series response (must be in (0,1))
#   link      : Character name of the link function (e.g., "logit")
#   parameters: List of model components and data (must include p, q, n, etc.)
#   ar        : Vector of AR lags
#   ma        : Vector of MA lags
#   X         : Optional exogenous variables (not used in current implementation)
#
# Output:
#   - Numeric vector of initial parameter estimates:
#       (alpha_start, varphi_start (AR), theta_start (MA), phi_start)

startValues <- function(y, link, parameters,
                        ar = NA, ma = NA, X = NA) {
  # from: 
  #       Beta Regression for Modelling Rates and Proportions
  #       Silvia Ferrari & Francisco Cribari-Neto
  #       p. 805
  
  # Parameters                                   -------------------------------
  
  
  # Parameters from outside this function
  p <- parameters$p
  q <- parameters$q
  n <- parameters$n
  
  p_length <- parameters$p_length
  q_length <- parameters$q_length
  
  ynew <- parameters$ynew
  
  scale <- parameters$scale
  
  # Parameters created locally
  m <- max(p, q, na.rm = TRUE)
  n <- length(y)
  
  linkfun <- parameters$linkfun
  linkinv <- parameters$linkinv
  mu.eta  <- parameters$mu.eta
  
  # P Matrix                                     -------------------------------
  
  P <- matrix(NA, nrow = n - m, ncol = p_length)
  for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]
  
  
  # Auxiliary Datasets                           -------------------------------
  
  x_inter <- matrix(1, nrow = n - m, ncol = 1)
  x_start <- cbind(x_inter, P)
  y_start <- linkfun(y[(m + 1):n])
  
  
  # Linear Regression for Initial Points         -------------------------------
  
  fit_start  <- lm.fit(x = x_start, y = y_start)
  
  mqo <- fit_start$coef
  
  alpha_start <- mqo[1]
  varphi_start <- mqo[-1]
  
  
  # Precision                                    -------------------------------
  
  k  <- length(mqo)
  n1 <- n - m
  
  y_hat_fit_start <- fitted(fit_start)
  mean_fit_start <- linkinv(y_hat_fit_start)
  
  linkfun_deriv <- 1 / mu.eta(eta = linkfun(mu = mean_fit_start))
  
  if (tolower(scale) == "predictor"){
    # Original link
    #er <- linkinv(y_start) - linkinv(y_hat_fit_start)
    er <- linkinv(residuals(fit_start))
  } else {
    # Predictor link
    er <- residuals(fit_start)
  }
  
  sigma2 <- sum(er^2) / ((n1 - k) * linkfun_deriv^2)
  
  phi_start_aux <- sum(mean_fit_start * (1 - mean_fit_start) / sigma2)
  phi_start <- phi_start_aux / n1
  
  # Initial value of \theta estimate
  theta_start <- rep(0, q_length)
  
  
  # Output                                       ------------------------------- 
  
  start_value <- c(alpha_start,
                   varphi_start, 
                   theta_start, 
                   phi_start)
  
  return(start_value)
}


# Log-likelihood - BARMA             -------------------------------------------

# Description:
#   Computes the negative log-likelihood for the Beta Autoregressive Moving Average (BARMA) model,
#   optionally allowing residuals to be computed on the predictor or response scale.
#
# Parameters:
#   z         : Numeric vector of current parameter estimates 
#               (alpha, varphi (AR), theta (MA), phi (precision))
#   parameters: List containing model and data settings including:
#                 - y        : original response (proportion)
#                 - ynew     : transformed response (for AR terms)
#                 - ar, ma   : AR and MA lags
#                 - scale    : either "predictor" or "original" for residual computation
#                 - linkfun, linkinv, mu.eta : link functions
#                 - n, m, p_length, q_length : numeric settings
#
# Output:
#   - Negative log-likelihood value (scalar)

loglik_arma <- function(z, parameters) {
  
  # Parameters          --------------------------------------------------------
  
  n        <- parameters$n
  p_length <- parameters$p_length
  q_length <- parameters$q_length
  
  m <- parameters$m
  
  y        <- parameters$y
  ynew     <- parameters$ynew
  
  ar       <- parameters$ar
  ma       <- parameters$ma
  
  scale <- parameters$scale
  
  alpha  <- z[1]
  varphi <- z[2:(p_length + 1)]
  theta  <- z[(p_length + 2):(p_length + q_length + 1)]
  phi    <- z[p_length + q_length + 2]
  
  linkfun <- parameters$linkfun
  linkinv <- parameters$linkinv
  mu.eta  <- parameters$mu.eta
  
  # Log-Likelihood      --------------------------------------------------------
  
  error <- rep(0, n)
  eta   <- rep(NA, n)
  mu    <- rep(NA, n)
  
  for (t in (m + 1):n) {
    eta[t]    <- alpha + varphi %*% ynew[t - ar] + theta %*% error[t - ma]
    mu[t]     <- linkinv(eta[t])
    error[t]  <- ifelse(tolower(scale) == "predictor",
                        ynew[t] - eta[t],
                        linkinv(ynew[t]) - mu[t])
  }
  
  mu_subset  <- mu[(m + 1):n]
  y_subset   <- y[(m + 1):n]
  
  ll_terms_arma <- dbeta(y_subset, mu_subset * phi, (1 - mu_subset) * phi,
                         log = TRUE)
  
  
  # Output              --------------------------------------------------------
  
  output <- (-1) * sum(ll_terms_arma)
  
  return(output)
}


# Score vector - BARMA               -------------------------------------------

# Description:
#   Computes the gradient (score vector) of the log-likelihood for the 
#   Beta Autoregressive Moving Average (BARMA) model. It supports both predictor-
#   and response-scale error formulations, and uses analytical derivatives
#   with respect to all parameters.
#
# Parameters:
#   z         : Numeric vector of current parameter estimates 
#               (alpha, varphi (AR), theta (MA), phi (precision))
#   parameters: List containing model components, such as:
#                 - y, ynew     : original and transformed response
#                 - ar, ma      : AR and MA lags
#                 - scale       : "predictor" or "original"
#                 - link, linkfun, linkinv, mu.eta : link function utilities
#                 - n, m, p_length, q_length       : dimensions
#
# Output:
#   - Score vector (gradient of the log-likelihood), negated for use with `optim()`

score_arma <- function(z, parameters) {
  
  # Parameters                        ------------------------------------------
  
  n        <- parameters$n
  m        <- parameters$m
  p_length <- parameters$p_length
  q_length <- parameters$q_length
  link     <- parameters$link
  
  y        <- parameters$y
  ynew     <- parameters$ynew
  
  ar       <- parameters$ar
  ma       <- parameters$ma
  
  link    <- parameters$link
  scale   <- parameters$scale
  
  alpha    <- z[1]
  varphi   <- z[2:(p_length + 1)]
  theta    <- z[(p_length + 2):(p_length + q_length + 1)]
  phi      <- z[p_length + q_length + 2]
  
  linkfun <- parameters$linkfun
  linkinv <- parameters$linkinv
  mu.eta  <- parameters$mu.eta
  
  
  # Regression Model for mu_t         ------------------------------------------
  
  error <- rep(0, n)
  eta   <- rep(NA, n)
  mu    <- rep(NA, n)
  
  for (t in (m + 1):n) {
    eta[t]    <- alpha + varphi %*% ynew[t - ar] + theta %*% error[t - ma]
    mu[t]     <- linkinv(eta[t])
    error[t]  <- ifelse(tolower(scale) == "predictor",
                        ynew[t] - eta[t],
                        linkinv(ynew[t]) - mu[t])
  }
  
  eta_subset <- eta[(m + 1):n]
  mu_subset  <- mu[(m + 1):n]
  y_subset   <- y[(m + 1):n]
  
  # Matrices R and P                  ------------------------------------------
  
  R <- matrix(nrow = n - m, ncol = q_length)
  for (t in 1:(n - m))  R[t, ] <- error[t + m - ma]
  
  P <- matrix(nrow = n - m, ncol = p_length)
  for (t in 1:(n - m)) P[t, ] <- ynew[t + m - ar]
  
  # Vector Score Initialization       ------------------------------------------
  
  mT     <- diag(mu.eta(eta = eta_subset))
  ystar  <- log(y_subset / (1 - y_subset))
  mustar <- digamma(mu_subset * phi) - digamma((1 - mu_subset) * phi)
  
  deta_dalpha  <- rep(0, n)
  deta_dvarphi <- matrix(0, nrow = n, ncol = p_length)
  deta_dtheta  <- matrix(0, nrow = n, ncol = q_length)
  deta         <- rep(NA, n)
  
  if (tolower(scale) == "original"){
    for (t in (m + 1):n) {
      mu_temp <- linkinv(eta_subset[t - ma])
      deta[t] <- inv_g_prime(mu_temp, link)
      
      deta_dalpha[t]  <- 1 - (deta[t] * theta %*% deta_dalpha[t - ma])
      deta_dvarphi[t, ] <- P[t - m, ] - (deta[t] * theta %*% deta_dvarphi[t - ma, ])
      deta_dtheta[t, ]  <- R[t - m, ] - (deta[t] * theta %*% deta_dtheta[t - ma, ])
    }
  } else {
    for (t in (m + 1):n) {
      deta_dalpha[t]  <- 1 - (theta %*% deta_dalpha[t - ma])
      deta_dvarphi[t, ] <- P[t - m, ] - (theta %*% deta_dvarphi[t - ma, ])
      deta_dtheta[t, ]  <- R[t - m, ] - (theta %*% deta_dtheta[t - ma, ])
    }
  }
  
  s  <- deta_dalpha[(m + 1):n]
  rP <- deta_dvarphi[(m + 1):n, ]
  rR <- deta_dtheta[(m + 1):n, ]
  
  # Vector Score                      ------------------------------------------
  
  ystar_mustar <- ystar - mustar
  mT_ystar_mustar <- crossprod(mT, ystar_mustar)
  
  U_alpha  <- phi * crossprod(s, mT_ystar_mustar)
  U_varphi <- phi * crossprod(rP, mT_ystar_mustar)
  U_theta  <- phi * crossprod(rR, mT_ystar_mustar)
  
  U_phi   <- sum(mu_subset * ystar_mustar + log(1 - y_subset) 
                 - digamma((1 - mu_subset) * phi) + digamma(phi))
  
  escore_vec <- c(U_alpha, U_varphi, U_theta, U_phi)
  
  # Output                            ------------------------------------------
  
  output <- (-1) * escore_vec
  
  return(output)
}


# Fisher information matrix - BARMA  -------------------------------------------

# Description:
#   Computes the observed Fisher information matrix for the BARMA model
#   (Beta Autoregressive Moving Average). It supports both the predictor-scale
#   and original-scale error formulation, and incorporates analytical derivatives.
#
# Parameters:
#   y         : Original time series (response values)
#   ar        : Vector of autoregressive lags
#   ma        : Vector of moving average lags
#   parameters: List with required elements:
#                 - y, ynew       : original and transformed series
#                 - n, m          : total length and initial burn-in
#                 - p_length, q_length : number of AR and MA parameters
#                 - link, linkfun, linkinv, mu.eta : link function details
#                 - scale         : model scale ("predictor" or "original")
#   alpha     : Intercept term
#   varphi    : Coefficients for autoregressive terms
#   theta     : Coefficients for moving average terms
#   phi       : Precision parameter
#   link      : Link function type (e.g., "logit", "probit", ...)
#
# Output:
#   A list with:
#     - fisher_info_mat: Fisher Information Matrix
#     - muhat          : Fitted mean values
#     - etahat         : Linear predictors
#     - errorhat       : Estimated residuals
#     - fitted         : Fitted ts object for plotting

inf_matrix_arma <- function(y, ar, ma, parameters,
                            alpha = 0, varphi = 0, theta = 0, 
                            phi = 0, link = link) {
  
  # Parameters                                              --------------------
  
  n        <- parameters$n
  m        <- parameters$m
  p_length <- parameters$p_length
  q_length <- parameters$q_length
  link     <- parameters$link
  
  y        <- parameters$y
  ynew     <- parameters$ynew
  
  link    <- parameters$link
  scale   <- parameters$scale
  
  linkfun <- parameters$linkfun
  linkinv <- parameters$linkinv
  mu.eta  <- parameters$mu.eta
  
  
  # Matrix P                                                --------------------
  
  P <- matrix(NA, nrow = n - m, ncol = p_length)
  for (i in 1:(n - m)) P[i, ] <- ynew[i + m - ar]
  
  
  # Estimation for error, mu, and eta                       --------------------
  
  errorhat <- rep(0, n)
  muhat <- rep(0, n)
  etahat <- rep(NA, n)
  
  for (t in (m + 1):n) {
    etahat[t] <- alpha + varphi %*% ynew[t - ar] + theta %*% errorhat[t - ma]
    muhat[t] <- linkinv(etahat[t])
    errorhat[t]  <- ifelse(tolower(scale) == "predictor",
                           ynew[t] - etahat[t],
                           linkinv(ynew[t]) - linkinv(etahat[t]))
  }
  
  etahat_subset <- etahat[(m + 1):n]
  muhat_subset <- muhat[(m + 1):n]
  y_subset <- y[(m + 1):n]
  
  
  # Matrix R                                                --------------------
  
  R <- matrix(nrow = n - m, ncol = q_length)
  for (i in 1:(n - m)) R[i, ] <- errorhat[i + m - ma]
  
  
  # Derivatives                                             --------------------
  
  deta_dalpha <- rep(0, n)
  deta.dvarphi <- matrix(0, ncol = p_length, nrow = n)
  deta_dtheta <- matrix(0, ncol = q_length, nrow = n)
  deta <- rep(0, n)
  
  if (tolower(scale) == "original"){
    for (i in (m + 1):n) {
      mu_temp <- linkinv(etahat_subset[i - ma])
      deta[i] <- inv_g_prime(mu_temp, link)
      
      deta_dalpha[i] <- 1 - (deta[i] * theta %*% deta_dalpha[i - ma])
      deta.dvarphi[i, ] <- P[(i - m), ] - (deta[i] * theta %*% deta.dvarphi[i - ma, ])
      deta_dtheta[i, ] <- R[(i - m), ] - (deta[i] * theta %*% deta_dtheta[i - ma, ])
    }
  } else {
    for (i in (m + 1):n) {
      deta_dalpha[i] <- 1 - (theta %*% deta_dalpha[i - ma])
      deta.dvarphi[i, ] <- P[(i - m), ] - (theta %*% deta.dvarphi[i - ma, ])
      deta_dtheta[i, ] <- R[(i - m), ] - (theta %*% deta_dtheta[i - ma, ])
    }
  }
  
  s <- deta_dalpha[(m + 1):n]
  rP <- deta.dvarphi[(m + 1):n, ]
  rR <- deta_dtheta[(m + 1):n, ]
  
  # Precompute some values                                  --------------------
  
  one_muhat <- 1 - muhat_subset
  
  psi1 <- trigamma(muhat_subset * phi)
  psi2 <- trigamma(one_muhat * phi)
  vc <- phi * (psi1 * muhat_subset - psi2 * one_muhat)
  D <- diag(psi1 * muhat_subset^2 + psi2 * one_muhat^2 - trigamma(phi))
  
  
  # Precompute vector mu_eta and matrix mT                  --------------------
  
  mu_eta <- mu.eta(eta = etahat_subset)
  mT <- diag(mu_eta)
  
  
  # Precompute W %*% s, W %*% rP, W %*% rR and mT %*% rR    --------------------
  
  W <- diag(phi * (psi1 + psi2)) %*% mT^2
  W_diag <- diag(W)
  
  W_s <- W_diag *  s
  W_rP <- W_diag * rP
  W_rR <- W_diag *  rR
  mT_vc <- mu_eta * vc
  
  
  # Compute the intermediate matrices using crossprod()     --------------------
  
  K_a_a <- phi * crossprod(s, W_s)
  K_p_a <- phi * crossprod(rP, W_s)
  K_t_a <- phi * crossprod(rR, W_s)
  
  K_p_p <- phi * crossprod(rP, W_rP)
  K_p_t <- phi * crossprod(rP, W_rR)
  K_t_t <- phi * crossprod(rR, W_rR)
  
  K_a_phi <- crossprod(s, mT_vc)
  K_p_phi <- crossprod(rP, mT_vc)
  K_t_phi <- crossprod(rR, mT_vc)
  K_phi_phi <- sum(diag(D))
  
  
  # Compute the remaining elements                          --------------------
  
  K_a_p <- t(K_p_a)
  K_t_p <- t(K_p_t)
  K_a_t <- t(K_t_a)
  K_phi_a <- K_a_phi
  K_phi_p <- t(K_p_phi)
  K_phi_t <- t(K_t_phi)
  
  
  # Construct the information matrix                        --------------------
  
  fisher_info_mat <- rbind(
    cbind(K_a_a, K_a_p, K_a_t, K_a_phi),
    cbind(K_p_a, K_p_p, K_p_t, K_p_phi),
    cbind(K_t_a, K_t_p, K_t_t, K_t_phi),
    cbind(K_phi_a, K_phi_p, K_phi_t, K_phi_phi)
  )
  
  if (any(is.na(ar)) == F) names_varphi <- paste("varphi", ar, sep = "")
  if (any(is.na(ma)) == F) names_theta  <- paste("theta", ma, sep = "")
  
  names_fisher_info_mat <- c("alpha", names_varphi, names_theta, "phi")
  colnames(fisher_info_mat) <- names_fisher_info_mat
  rownames(fisher_info_mat) <- names_fisher_info_mat
  
  
  # Output                                                  --------------------
  
  output_list <- list()
  
  output_list$fisher_info_mat <- fisher_info_mat
  
  output_list$muhat <- muhat
  output_list$etahat <- etahat
  output_list$errorhat <- errorhat
  
  # fitted values 
  fitted_values <- ts(c(rep(NA, m), muhat[(m + 1):n]), 
                      start = start(y), 
                      frequency = frequency(y))
  
  output_list$fitted <- fitted_values
  
  return(output_list)
}


# Last touches / Model output        -------------------------------------------

# Description:
#   This function generates the final model output, including parameter estimates,
#   standard errors, confidence metrics, and information criteria, based on the
#   fitted BARMA model and its Fisher information matrix.
#
# Parameters:
#   result           : List containing the model output from estimation
#   fisher_info_mat  : Observed Fisher Information Matrix
#   h                : Forecast horizon
#
# Global dependencies:
#   - parameters     : List with n, m, p_length, q_length, link
#   - start_values   : Vector of initial estimates
#   - y_prev         : Vector with extended series (e.g., including forecasts)
#   - X              : Optional design matrix for exogenous covariates
#   - beta           : Optional coefficient vector for covariates
#
# Output:
#   The input 'result' list is updated with:
#     - vcov                : Inverse Fisher Information matrix
#     - forecast            : Forecasted values (ts)
#     - model               : Table of estimates, std. errors, z-stats, p-values
#     - aic, bic, hq        : Information criteria

final_output <- function(parameters, result, X, start_values, fisher_info_mat, h){
  
  # Extract parameters from global environment
  n        <- parameters$n
  m        <- parameters$m
  
  p_length <- parameters$p_length
  q_length <- parameters$q_length
  
  link     <- parameters$link
  
  y_prev   <- parameters$y_prev
  
  # Inverse of Fisher information matrix
  vcov <- try(solve(fisher_info_mat, tol = 1e-20), silent = TRUE)
  
  # Check if the inverse matrix is possible
  if (!(typeof(vcov) == "double")) {
    warning("FISHER'S INFORMATION MATRIX IS NOT INVERTIBLE! ")
  } else {
    
    result$inv_inf_matrix <- 1 # Reamostragem
    #result$inv_inf_matrix <- 0
    
    # output
    result$start_values <- start_values
    result$fisher_info_mat <- fisher_info_mat
    result$forecast <- y_prev[(n + 1):(n + h)]
    
    # output inverse of Fisher information matrix
    result$vcov <- vcov
    
    # Model presentation ----
    stderror <- sqrt(diag(vcov))
    z_stderror <- stderror
    
    z_zstat <- abs(result$coef / stderror)
    z_pvalues <- 2 * (1 - pnorm(z_zstat))
    
    model_presentation <- cbind(round(result$coef, 4),
                                round(z_stderror, 4),
                                round(z_zstat, 4),
                                round(z_pvalues, 4))
    
    colnames(model_presentation) <- c("Estimate",
                                      "Std. Error",
                                      "z value",
                                      "Pr(>|z|)")
    
    result$model <- model_presentation
    
    result$link <- link
    
    # information criteria
    if (any(is.na(X) == FALSE)) {
      
      aux_info1 <- -2 * result$loglik
      aux_info2 <- p_length + q_length + 2 + length(beta)
      log_n <- log(n)
      
      result$aic <- aux_info1 + 2 * aux_info2
      result$bic <- aux_info1 + log_n * aux_info2
      result$hq  <- aux_info1 + log(log_n) * 2 * aux_info2
      
    } else {
      
      aux_info1 <- -2 * result$loglik
      aux_info2 <- p_length + q_length + 2
      log_n <- log(n)
      
      result$aic <- aux_info1 + 2 * aux_info2
      result$bic <- aux_info1 + log_n * aux_info2
      result$hq  <- aux_info1 + log(log_n) * 2 * aux_info2
      
    }
  }
  
  return(result)
}
