#' ============================================================================
#' Code for performing simulation data from \beta ARMA model
#' ============================================================================
#'  
#' -------------------------------------------------------------------------- #
#' AUTHOR
#' -------------------------------------------------------------------------- #
#' Based on the original code by Fabio M. Bayer (bayer@ufsm.br) 
#' on 2015-10-15
#' 
#' Modified and improved by Everton da Costa (everton.ecosta@ufpe.br)
#' on 2022-02-03
#' 
#' VERSION: 2.00
#' IDENTIFIER: simuBarma 
#' UPDATE: 2023-02-03
#' 
#' 
#'--------------------------------------------------------------------------- #
#' DESCRIPTION:
#'--------------------------------------------------------------------------- #
# This code implements a beta ARMA analysis, building upon the original code
# developed by Fabio M Bayer. The modifications made by Everton da Costa
# include:
# - Rewriting all the code to follow the tidyverse style guide
# Usage:
#'--------------------------------------------------------------------------- #
#' INPUT
#'--------------------------------------------------------------------------- #
#' @param n The length of the time series to be simulated
#' @param varphi A vector of autoregressive (AR) parameters
#' @param theta A vector of moving average (MA) parameters
#' @param alpha The intercept term
#' @param phi The precision parameter of the \beta ARMA distribution
#' @param freq The frequency of the time series (e.g., 12 for monthly data)
#' @param link The link function to be used ("logit", "loglog", "cloglog")
#'
#'
#'--------------------------------------------------------------------------- #
#' OUTPUT
#'--------------------------------------------------------------------------- #
#' @return A time series object containing the simulated \beta ARMA process
#' 
# ===========================================================================

simuBarma <- function(n,
                      varphi = NA, theta = NA,
                      alpha = 0.0, phi = 20, 
                      freq = 12, link = "logit") {
  
  ar <- NA
  ma <- NA
  
  if (any(is.na(varphi) == F)) ar <- 1:length(varphi)
  if (any(is.na(theta) == F)) ma <- 1:length(theta)
  
  # Link functions
  # The link function structure provides the necessary components
  # to perform the \beta ARMA model simulation with the specified link function.
  # Available link functions: "logit", "loglog", "cloglog", "probit"
  # ---------------------------------------------------------------------------
  # Link functions
  # ---------------------------------------------------------------------------
  
  # Check if the link argument is a character string or an expression
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
    if (linktemp == "link") {
      linktemp <- eval(link)
    }
  }
  
  # Set up the link function based on the provided link argument
  if (linktemp == "logit") {
    
    # Use make.link function for logit link (efficient C implementation)
    stats <- make.link("logit")
    
  } else if (linktemp == "probit") {
    
    # Use make.link function for probit link
    stats <- make.link("probit")
    
  } else if (linktemp == "cloglog") {
    
    # Use make.link function for cloglog link
    stats <- make.link("cloglog")
    
  } else if (linktemp == "loglog") {
    
    # Manually define linkfun, linkinv, and mu.eta for loglog link
    stats <- list()
    
    stats$linkfun <- function(mu) -log(-log(mu))
    stats$linkinv <- function(eta) exp(-exp(-eta))
    stats$mu.eta <- function(eta) exp(-exp(-eta)) * exp(-eta)
    
  } else {
    # If the provided link is not supported, raise an error
    stop(
      paste0(
        "The '", linktemp, "' link function is not available.", "\n",
        "Available options are: ", "logit, loglog, cloglog, and probit."
      )
    )
  }
  
  # Create a list with the link function details
  link1 <- structure(list(
    link = linktemp,
    linkfun = stats$linkfun,
    linkinv = stats$linkinv,
    mu.eta = stats$mu.eta
  ))
  
  # Assign link function components to separate variables
  linkfun <- link1$linkfun
  linkinv <- link1$linkinv
  mu.eta <- link1$mu.eta
  
  # ===========================================================================
  # ARMA model
  # ===========================================================================
  if (any(is.na(varphi) == F) && any(is.na(theta) == F)) {
    # This section simulates a \beta ARMA(p ,q) process
    # where both AR and MA orders are specified
    
    # Compute the AR and MA orders
    p <- max(ar)
    q <- max(ma)
    m <- max(p, q)
    
    ynew <- rep(alpha, (n + m))
    mu <- linkinv(ynew)
    
    # E(error) = 0
    error <- rep(0, n + m)  
    eta <- y <- rep(NA, n + m)
    
    for (i in (m + 1):(n + m)) {
      
      eta[i] <- alpha + varphi %*% ynew[i - ar] + theta %*% error[i - ma]
      mu[i] <- linkinv(eta[i])
      y[i] <- rbeta(1, mu[i] * phi, (1 - mu[i]) * phi)
      ynew[i] <- linkfun(y[i])
      error[i] <- ynew[i] - eta[i]
      
    }
    
    return(ts(y[(m + 1):(n + m)], frequency = freq))
  }
  
  
  # ===========================================================================
  # AR model
  # ===========================================================================
  if (any(is.na(varphi) == F) && any(is.na(theta) == T)) {
    # This section simulates a \beta AR(p) process
    # where only the AR order is specified
    
    # Compute the AR order
    p <- max(ar)
    m <- p
    
    ynew <- rep(alpha, (n + m))
    mu <- linkinv(ynew)
    
    eta <- y <- rep(NA, n + m)
    
    for (i in (m + 1):(n + m)) {
      
      eta[i] <- alpha + varphi %*% ynew[i - ar]
      mu[i] <- linkinv(eta[i])
      y[i] <- rbeta(1, mu[i] * phi, (1 - mu[i]) * phi)
      ynew[i] <- linkfun(y[i])
      
    }
    
    return(ts(y[(m + 1):(n + m)], frequency = freq))
  }
  
  
  # =========================================================================== 
  # MA model
  # ===========================================================================
  if (any(is.na(varphi) == T) && any(is.na(theta) == F)) {
    # This section simulates a \beta MA(q) process
    # where only the MA order is specified
    
    # Compute the MA order
    q <- max(ma)
    m <- q
    
    ynew <- rep(alpha, (n + m))
    mu <- linkinv(ynew)
    
    # E(error)=0
    eta <- y <- error <- rep(0, n + m)
    
    for (i in (m + 1):(n + m)) {
      
      eta[i] <- alpha + theta %*% error[i - ma]
      mu[i] <- linkinv(eta[i])
      y[i] <- rbeta(1, mu[i] * phi, (1 - mu[i]) * phi)
      ynew[i] <- linkfun(y[i])
      error[i] <- ynew[i] - eta[i]
      
    }
    
    return(ts(y[(m + 1):(n + m)], frequency = freq))
  }
}