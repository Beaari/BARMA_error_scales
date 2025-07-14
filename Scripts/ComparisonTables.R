# =============================================================================
# Comparison of estimated parameters for BARMA with the error in the predictor
# scale against the error on the original scale
#
# Date: 2025/07/13
# 
# =============================================================================

# Load required packages                ----------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  zoo,      # For working with yearmon dates
  dplyr,    # For piping the data
  ggplot2,  # For plotting in a tidy manner
  Metrics,  # For calculating RMSE
  moments,  # For calculating moments, skewness and kurtosis
  forecast  # For time series forecasting
)


# Load custom functions                 ----------------------------------------

options(scipen = 5)

source("Scripts/Functions/SimulationBarma.R")    # Function to simulate BARMA data

source("Scripts/Functions/BARMA.R")              # BARMA function
source("Scripts/Functions/AuxiliaryFunctions.R") # Auxiliary functions to BARMA model

get_simulated_tab <- function(settings_par, seed, h){
  
  # Settings of data generation           --------------------------------------
  
  sample_size <- settings_par[[1]]   # Set the number of observations to simulate
  alpha_true  <- settings_par[[2]]   # Set the intercept for the BARMA model
  varphi_true <- settings_par[[3]]   # Set the autoregressive parameter
  theta_true  <- settings_par[[4]]   # Set the moving average parameter
  phi_true    <- settings_par[[5]]   # Set the precision parameter
  link        <- settings_par[[6]]   # Specify the link function for the model
  
  
  # Simulate Data                         --------------------------------------
  
  set.seed(seed)
  
  # Generate data using the specified BARMA parameters
  y_orig <- simuBarma(n = sample_size, 
                 varphi = varphi_true, 
                 theta = theta_true, 
                 alpha = alpha_true, 
                 phi = phi_true, 
                 link = link,
                 scale = "original")
  
  y_pred <- simuBarma(n = sample_size, 
                      varphi = varphi_true, 
                      theta = theta_true, 
                      alpha = alpha_true, 
                      phi = phi_true, 
                      link = link,
                      scale = "predictor")
  
  # Fit BARMA Model                       --------------------------------------
  
  ar_vec <- settings_par[[7]]        # Define autoregressive order
  ma_vec <- settings_par[[8]]        # Define moving average order
  
  y_window_orig <- window(y_orig, end = time(y_orig)[length(y_orig) - h])
  y_window_pred <- window(y_pred, end = time(y_pred)[length(y_pred) - h])
  
  # Fit the BARMA with error on original scale with chosen link
  fit_BARMA_orig <- barma(y = y_window_orig, 
                          ar = ar_vec,
                          ma = ma_vec,
                          link = link,
                          scale = "original",
                          h = h,
                          X = NA, 
                          X_hat = NA)
  
  # Fit the BARMA with error on predictor with chosen link
  fit_BARMA_pred <- barma(y = y_window_pred, 
                          ar = ar_vec,
                          ma = ma_vec,
                          link = link,
                          scale = "predictor",
                          h = h,
                          X = NA, 
                          X_hat = NA)
  
  
  # Metric                                --------------------------------------
  
  y_obs_orig <- window(y_orig, start = time(y_orig)[length(y_orig) - h + 1])
  y_obs_pred <- window(y_pred, start = time(y_pred)[length(y_pred) - h + 1])
  
  RMSE_orig <- rmse(y_obs_orig, fit_BARMA_orig$forecast) # Original scale
  RMSE_pred <- rmse(y_obs_pred, fit_BARMA_pred$forecast) # Predictor scale
  
  
  # Output                                --------------------------------------
  
  return(list(y_orig = y_orig,
              y_pred = y_pred,
              fit_BARMA_pred = fit_BARMA_pred, 
              fit_BARMA_orig = fit_BARMA_orig,
              RMSE_pred = RMSE_pred,
              RMSE_orig = RMSE_orig))
}


# Simulate Data and fit model           ----------------------------------------

# Settings of data generation:
# In order: n, alpha, varphi, theta, phi, link, ar, ma
settings_par <- list(240, 0, 0.4, -0.7, 20, "logit", 1, 1)


objs <- get_simulated_tab(settings_par, seed = 16, h = 6)

y <- objs$y_orig

# Data behavior through time-series, ACF, and PACF plots
y %>% 
  ggtsdisplay(
    theme = theme_light()+
      theme(panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank()))


# Example: Fit BARMA(1, 1) Model        ----------------------------------------

fit_BARMA_pred <- objs$fit_BARMA_pred # Predictor scale

fit_BARMA_orig <- objs$fit_BARMA_orig # Original scale

# Convergence successful
fit_BARMA_pred$conv == 0
fit_BARMA_orig$conv == 0

# Estimates of the model
fit_BARMA_pred$coef
fit_BARMA_orig$coef

# RMSE
objs$RMSE_pred
objs$RMSE_orig

# Comparison table                      ----------------------------------------

list_settings_par <- list(
  settings_par = list(240, 0, 0.29, -0.25, 15, "logit", 1, 1),
  settings_par = list(240, 0, 0.29, 0.25, 15, "logit", 1, 1)
)

df_results <- lapply(list_settings_par, function(x){
  objs <- get_simulated_tab(x, seed = 16, h = 6)
  
  data.frame(
    true_parameters = c(unlist(x)[2:5],"AIC","BIC","RMSE"),
    orig_parameters = c(objs$fit_BARMA_orig$coef,
                        objs$fit_BARMA_orig$aic,
                        objs$fit_BARMA_orig$bic,
                        objs$RMSE_orig),
    pred_parameters = c(objs$fit_BARMA_pred$coef,
                        objs$fit_BARMA_pred$aic,
                        objs$fit_BARMA_pred$bic,
                        objs$RMSE_pred)
  ) %>%
    mutate(across(c("orig_parameters","pred_parameters"), ~round(.x, 4)))
})

df_results[[1]]
df_results[[2]]
