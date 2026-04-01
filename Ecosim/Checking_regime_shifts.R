# Load required packages
library(forecast)
library(nlme)
#install.packages('rshift', repos = c('https://alexhroom.r-universe.dev', 'https://cloud.r-project.org'), type='binary' )

library(rshift)
library(dplyr)
library(lubridate)
library(gridGraphics)
library(ggplot2)
library(MASS)

## function to get regime shifts using package 'rshift'
get_regime_shifts <- function(df,l){## Calculate Anomalies in Time Series Data with rshift (STARS)
  plot_list <- list()
  # ecoind_mc
  for (i in colnames(df[-1])) {
    # Step 1: Detrend using OLS
    #model_ols <- lm(df[,i] ~  df$date)
    model_ols <- rlm(df[,i] ~ df$date, data = df, method = "M", psi = psi.huber, k2 = 3) #, k2 = 3 adjusting weight of outliers.

    detrended_data <- residuals(model_ols)
    # Plot detrended data
    plot(df$date, detrended_data, type = "l", col = "red", main = "Detrended Data")

    # Step 2: Account for autocorrelation using ARIMA
    arima_model <- auto.arima(detrended_data)
    anomalies <- residuals(arima_model)

    # Step 3: Visualize anomalies
    plot(anomalies, type = "l", col = "blue", main = "Anomalies Over Time")


    # STARS (Sequential T-test Analysis of Regime Shifts)
    anomaly_data <- cbind(df[,c("date",i)],anomalies)

    RSI_data <- Rodionov(data=anomaly_data, col="anomalies" , time= "date", l=l, prob = 0.05,
                         startrow = 1, merge = TRUE)

    #ecoind_NA_annual$time <- 1:33
    #Lanzante(data=ecoind_NA_annual, col="Biomass" , time= "time", p = 0.05, merge = FALSE)

    RSI_graph(RSI_data, "anomalies", "date", "RSI", mean_lines = TRUE)
    plot_list[[i]] <- recordPlot()


  }
  return(plot_list)
}

get_regime_shifts_raw <- function(df,l){## Calculate Anomalies in Time Series Data with rshift (STARS)
  plot_list_raw <- list()
  # ecoind_mc
  for (i in colnames(df[-1])) {
    # STARS (Sequential T-test Analysis of Regime Shifts)
    RSI_data <- Rodionov(data= df, col= i , time= "date", l=l, prob = 0.05,
                         startrow = 1, merge = TRUE)

    #ecoind_NA_annual$time <- 1:33
    #Lanzante(data=ecoind_NA_annual, col="Biomass" , time= "time", p = 0.05, merge = FALSE)

    RSI_graph(RSI_data, i, "date", "RSI", mean_lines = TRUE)
    plot_list_raw[[i]] <- recordPlot()


  }
  return(plot_list_raw)
}

## prepare eco indicators for regime shift check
# take annual values
ecoind_sim$Time <- ecoind_mc_grouped$Time
ecoind_regimeshift <- as.data.frame(ecoind_sim %>%
  mutate(date = year(Time)) %>%
  group_by(date) %>%
  summarize( med_total_B = mean(Total.B),
             med_total_C = mean(Total.C),
             med_TL_catch = mean(TL.catch),
             med_TL_community = mean(TL.community)))



## get regime shift for indicators ran with MC
plots_mc <- get_regime_shifts(df = ecoind_regimeshift,l = 10)

plots_mc$med_total_B

plots_mc$med_total_C

plots_mc$med_TL_catch

plots_mc$med_TL_community

## get regime shift for indicators ran with MC
plots_na <- get_regime_shifts(df = ecoind_NA_annual, l = 10)

plots_na$TST

plots_na$IFO

plots_na$FCI

plots_na$OC

plots_na$AC



# raw data
## get regime shift for indicators ran with MC
plots_mc_raw <- get_regime_shifts_raw(df = ecoind_regimeshift,l = 10)

plots_mc_raw$med_total_B

plots_mc_raw$med_total_C

plots_mc_raw$med_TL_catch

plots_mc_raw$med_TL_community

## get regime shift for indicators ran with MC
plots_na_raw <- get_regime_shifts_raw(df = ecoind_NA_annual, l = 10)

plots_na_raw$TST

plots_na_raw$IFO

plots_na_raw$FCI

plots_na_raw$OC

plots_na_raw$AC

