# =============================================================================
# Example usage of the BARMA model using ONS hidrological data
#
# Date: 2025/07/13
# 
# =============================================================================

## Packages                      -----------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  zoo,      # For working with yearmon dates
  urca,     # For testing stationarity
  dplyr,    # For piping the data
  aws.s3,   # For accessing ONS data
  ggplot2,  # For plotting in a tidy manner
  Metrics,  # For calculating RMSE
  moments,  # For calculating moments, skewness and kurtosis
  forecast, # For time series forecasting
  openxlsx, # For reading .xlsx files
  lubridate # For working with dates
)


## Functions                     -----------------------------------------------

source("Scripts/Functions/BARMA.R")              # BARMA function
source("Scripts/Functions/AuxiliaryFunctions.R") # Auxiliary functions to BARMA model

prediction_plot <- function(fit_B,n,h){
  ons_volume_df %>% 
    select(date, mean_volutil) %>% 
    mutate(pred = "N") %>%
    bind_rows(
      data.frame(
        date = ons_volume_df$date[n:dim(ons_volume_df)[1]],
        mean_volutil = fit_B$forecast,
        pred = rep("Y",h))
    ) %>% 
    filter(date >= "2020-01-01") %>% 
    ggplot()+
    geom_line(aes(x=date, y = mean_volutil, color = pred), size = 1)+
    theme_light()
}


## Loading ONS data from AWS     -----------------------------------------------

# Commenting this section out, this takes awhile to run, if needed to
# replicate just uncomment it, it takes about 2 ~ 3 minutes to run all
# Below we run a .RDS file with all the collected data for quick use

#bucket <- "ons-aws-prod-opendata" # AWS bucket for ONS data
#
#Sys.setenv("AWS_DEFAULT_REGION" = "us-west-2") # Set correct region
#
## Gather bucket names for datasets we'll use
#objs <- get_bucket(bucket = bucket, 
#                   prefix = "dataset/dados_hidrologicos_di/DADOS_HIDROLOGICOS",
#                   max = Inf)
#
## Keys with the path in AWS
#keys <- data.frame(key = sapply(objs, function(x) x[["Key"]])) %>% 
#  filter(grepl(".xlsx",key)) %>% 
#  pull(key)
#
#ons_volume_list <- lapply(keys, function(x){
#  # Reading ONS data from AWS bucket
#  ons_volume_s3 <- s3read_using(
#    FUN = read.xlsx,
#    object = x,
#    bucket = bucket
#  )
#  
#  ons_volume_s3 <- ons_volume_s3 %>% 
#    select(din_instante,nom_reservatorio,val_volumeutilcon)
#  
#  return(ons_volume_s3)
#})
#
#ons_volume_df <- bind_rows(ons_volume_list) # Binding the list

# Saving ons data for quick access
#saveRDS(ons_volume_df,"Data/ons_volume_df.RDS")

ons_volume_df <- readRDS("Data/ons_volume_df.RDS")


## Cleaning ONS data             -----------------------------------------------

# Transforming date and filtering for C.BRANCO-1 reservoir
ons_volume_df <- ons_volume_df %>% 
  as_tibble() %>% 
  mutate(date = as.numeric(din_instante),
         date = as.Date(date, origin = "1899-12-30"),
         val_volumeutilcon = val_volumeutilcon/100) %>%
  filter(date >= "2016-01-01") %>% 
  filter(nom_reservatorio == "C.BRANCO-1") %>%
  select(nom_reservatorio,date,val_volumeutilcon)

# Filtering 2016 onward and calculating mean per month and year
ons_volume_df <- ons_volume_df %>% 
  mutate(year = year(date),
         month = month(date)) %>% 
  group_by(year,month) %>% 
  summarise(
    date = min(date),
    mean_volutil = mean(val_volumeutilcon, na.rm = T))


## Time-Series Behavior          -----------------------------------------------

# Defining time-series
y_ts <- ts(ons_volume_df$mean_volutil, freq = 12)

# Plotting time-series correlograms
ggtsdisplay(ts(ons_volume_df$mean_volutil, freq = 12), 
            theme = theme_light()+theme(panel.grid.minor = element_blank(),
                                        panel.grid.major.x = element_blank(),))

summary(y_ts) # Summary stats
mean(y_ts)    # Mean of the time series
sd(y_ts)      # Standard Deviation of the series

# We failed to reject the null hypothesis of stationarity
# This series seems to be stationary
ur.fit = ur.df(y_ts, type="drift", lags=12, selectlags="AIC")
summary(ur.fit) 


## Fitting BARMA to ONS data     -----------------------------------------------

h <- 6
num <- length(y_ts) - h

# Fitting BARMA with error on the original scale and linking 
# function 'logit'
fit_BARMA_orig <- barma(y = y_ts, 
                    ar = 1,
                    ma = 1,
                    link = "logit",
                    scale = "original",
                    h = h,
                    X = NA, 
                    X_hat = NA)

# Fitting BARMA with error on the predictor scale and linking 
# function 'logit'
fit_BARMA_pred <- barma(y = y_ts, 
                    ar = 1,
                    ma = 1,
                    link = "logit",
                    scale = "predictor",
                    h = h,
                    X = NA, 
                    X_hat = NA)

fit_BARMA_orig$coef  # Parameters estimations / Original Scale
fit_BARMA_pred$coef  # Parameters estimations / Predictor Scale


## Forecasting                   -----------------------------------------------

prediction_plot(fit_BARMA_orig,num+1,h) # Predictions with original scale
prediction_plot(fit_BARMA_pred,num+1,h) # Predictions with predictor scale

# Creating a dataset with the output
df_preds <- data.frame(
  date = ons_volume_df$date[(num+1):dim(ons_volume_df)[1]],
  observed = y_ts[(num+1):length(y_ts)],
  prev_original = fit_BARMA_orig$forecast,
  prev_predictor = fit_BARMA_pred$forecast
  ) %>% 
  mutate(across(c("observed",starts_with("prev")), ~round(.x,4)))

df_preds %>% 
  mutate(across("observed", ~as.character(.x))) %>% 
  bind_rows(
    data.frame(
      date = c(NA,NA,NA),
      observed = c("RSME","AIC","BIC"),
      prev_original = c(rmse(df_preds$observed,df_preds$prev_original),
                        fit_BARMA_orig$aic,fit_BARMA_orig$bic),
      prev_predictor = c(rmse(df_preds$observed,df_preds$prev_predictor),
                         fit_BARMA_pred$aic,fit_BARMA_pred$bic)
    ) %>% 
      mutate(across(starts_with("prev"), ~round(.x,4)))
  )
     