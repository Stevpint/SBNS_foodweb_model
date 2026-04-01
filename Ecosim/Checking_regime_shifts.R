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
## function RSI_graph adjusted
RSI_graph_SP <- function (data, col, time, rsi, mean_lines = FALSE){
  p1 <- ggplot(data) + geom_col(aes(x = .data[[time]], y = .data[[col]]))+
      theme_bw()+
      theme(
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 90, hjust = 1, colour = "black", size = 14),
        axis.text.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 16,),
        axis.title.y = element_text(colour = "black", size = 16,),
        strip.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
      )+
    labs(x = 'Year', y = "Anomalies")
  if (mean_lines) {
    means <- rshift::regime_means(data, col, rsi)
    p1 <- p1 + geom_line(aes(x = .data[[time]], y = means),
                         color = "red")
  }
  p2 <- ggplot(data) + geom_col(aes(x = .data[[time]], y = .data[[rsi]]))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(2,
                                                               1)))
  vplayout <- function(x, y) grid::viewport(layout.pos.row = x,
                                            layout.pos.col = y)
  print(p1, vp = vplayout(1, 1))
  print(p2, vp = vplayout(2, 1))
  return(p1)
}
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



    plot_list[[i]] <- RSI_graph_SP(RSI_data, "anomalies", "date", "RSI", mean_lines = TRUE)


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

    plot_list_raw[[i]] <- RSI_graph_SP(RSI_data, i, "date", "RSI", mean_lines = TRUE)


  }
  return(plot_list_raw)
}

## prepare eco indicators for regime shift check
# take annual values
ecoind_sim$Time <- ecoind_mc_grouped$Time
ecoind_regimeshift <- as.data.frame(ecoind_sim %>%
  mutate(date = year(Time)) %>%
  group_by(date) %>%
  summarize( total_B = mean(Total.B),
             total_C = mean(Total.C),
             TL_catch = mean(TL.catch),
             TL_community = mean(TL.community)))



## get regime shift for indicators ran with MC
plots_mc <- get_regime_shifts(df = ecoind_regimeshift,l = 10)

plots_mc$total_B <- plots_mc$total_B + ggtitle("TBco") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plots_mc$total_C <- plots_mc$total_C + ggtitle("TC") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plots_mc$TL_catch <- plots_mc$TL_catch + ggtitle("TLc") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plots_mc$TL_community <- plots_mc$TL_community + ggtitle("mTLco") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

## get regime shift for indicators ran with MC
plots_na <- get_regime_shifts(df = ecoind_NA_annual, l = 10)

plots_na$TST <- plots_na$TST + ggtitle("TST") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plots_na$IFO <- plots_na$IFO + ggtitle("IFO") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plots_na$FCI <- plots_na$FCI + ggtitle("FCI") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plots_na$OC <- plots_na$OC + ggtitle("O/C") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plots_na$AC <- plots_na$AC + ggtitle("A/C") +
  theme(
    plot.title = element_text(
      hjust = 0.5,  # Center the title
      size = 16,    # Set font size
      face = "bold" # Make it bold
    )
  )

plot_grid(plots_mc$total_B,
          plots_mc$total_C,
          plots_mc$TL_catch,
          plots_mc$TL_community,
          plots_na$TST,
          plots_na$IFO,
          plots_na$FCI,
          plots_na$OC,
          plots_na$AC, labels = NA, ncol = 3, align = 'v') +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


# raw data (not used)
## get regime shift for indicators ran with MC
plots_mc_raw <- get_regime_shifts_raw(df = ecoind_regimeshift,l = 10)

plots_mc_raw$total_B

plots_mc_raw$total_C

plots_mc_raw$TL_catch

plots_mc_raw$TL_community

## get regime shift for indicators ran with MC
plots_na_raw <- get_regime_shifts_raw(df = ecoind_NA_annual, l = 10)

plots_na_raw$TST

plots_na_raw$IFO

plots_na_raw$FCI

plots_na_raw$OC

plots_na_raw$AC

