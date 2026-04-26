# =====================================================
# PROJECT: Durham Temperature Forecasting
# AUTHOR: David A. Olaechea Dongo
# DESCRIPTION:
# Forecast daily mean temperature in Durham using 
# 1901–2019 data and evaluate against 2020 observations.
# =====================================================

#set working directory
setwd("C:/R_Projects/DAES/Mini-Project")

# Load libraries
library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(zoo)
library(Metrics)
library(gridExtra)

# Load the daily data 1901-2019
DOCtemp <- read.csv("durhamtemp_1901_2019.csv")

# Convert date
DOCtemp$Date <- dmy(DOCtemp$Date)

# check conversion
sum(is.na(DOCtemp$Date))
str(DOCtemp$Date)

# order by date
DOCtemp <- DOCtemp[order(DOCtemp$Date), ]

# Missing dates
full_dates <- seq.Date(min(DOCtemp$Date), max(DOCtemp$Date), by = "day")
missing_dates <- setdiff(full_dates, DOCtemp$Date)
n_missing_dates <- length(missing_dates)

n_missing_dates
head(missing_dates)

# Duplicate dates
n_duplicate_dates <- sum(duplicated(DOCtemp$Date))
n_duplicate_dates

# Missing values
na_summary <- colSums(is.na(DOCtemp))
na_summary

# Raw daily temperature with smooth trend
ggplot(DOCtemp, aes(x = Date, y = Av.temp)) +
  geom_line(colour = "grey40", linewidth = 0.2) +
  geom_smooth(method = "loess", se = FALSE, colour = "firebrick", linewidth = 1) +
  xlab("Date") +
  ylab("Daily average temperature (°C)") +
  labs(title="Average Temperature Time Series (1901 - 2019)") +
  theme_minimal()

# Month factor for plotting
DOCtemp$Month_f <- factor(
  DOCtemp$Month,
  levels = 1:12,
  labels = c("Jan","Feb","Mar","Apr","May","Jun",
             "Jul","Aug","Sep","Oct","Nov","Dec")
)

# Monthly boxplot
ggplot(DOCtemp, aes(x = Month_f, y = Av.temp)) +
  geom_boxplot() +
  xlab("Month") +
  ylab("Daily average temperature (°C)") +
  labs(title="Monthly Boxplots for Daily Average Temperature") +
  theme_minimal()

# 30-day moving average
rolling_30 <- DOCtemp %>%
  mutate(
    roll30 = rollmean(Av.temp, k = 30, fill = NA, align = "center")
  )

c1 <- ggplot(rolling_30, aes(x = Date, y = roll30)) +
  geom_line(colour = "black", linewidth = 0.5) +
  xlab("Date") +
  ylab("30-day moving average of daily temperature (°C)") +
  labs(title = "30-Day Moving Average") +
  theme_minimal()

# 365-day moving average
rolling_365 <- DOCtemp %>%
  mutate(
    roll365 = rollmean(Av.temp, k = 365, fill = NA, align = "center")
  )

c2 <- ggplot(rolling_365, aes(x = Date, y = roll365)) +
  geom_line(colour = "black", linewidth = 0.5) +
  xlab("Date") +
  ylab("365-day moving average of daily temperature (°C)") +
  labs(title = "365-Day Moving Average") +
  theme_minimal()

grid.arrange(c1, c2, ncol = 2)

# Stationarity and linear detrending

# ADF test
adf_result <- adf.test(DOCtemp$Av.temp)
adf_result

# Create numeric time index
DOCtemp$time_index <- 1:nrow(DOCtemp)

# Linear trend
model_linear <- lm(Av.temp ~ time_index, data = DOCtemp)
DOCtemp$linear_residual <- resid(model_linear)

# Plot fitted linear trend
ggplot(DOCtemp, aes(x = Date, y = Av.temp)) +
  geom_line(colour = "black", linewidth = 0.5) +
  stat_smooth(method = "lm", aes(colour = "Linear"), se = FALSE, linewidth = 0.8) +
  xlab("Date") +
  ylab("Daily average temperature (°C)") +
  labs(colour = "Trend model") +
  theme_minimal()

# Variance comparison
detrend <- data.frame(
  original_var = var(DOCtemp$Av.temp, na.rm = TRUE),
  linear_var = var(DOCtemp$linear_residual, na.rm = TRUE)
)

detrend

# Daily time series
temp_ts <- ts(DOCtemp$Av.temp, frequency = 365)

# Model 1: Seasonal naive
model_snaive <- snaive(temp_ts, h = 366)

summary(model_snaive)

autoplot(model_snaive) +
  xlab("Time") +
  ylab("Daily average temperature (°C)") +
  theme_minimal()

#deseasoning
DOCtemp$month_day <- format(DOCtemp$Date, "%m-%d")

seasonal_summary <- DOCtemp %>%
  group_by(month_day) %>%
  summarise(
    median_temp = median(Av.temp, na.rm = TRUE),
    upper = quantile(Av.temp, 0.75, na.rm = TRUE),
    lower = quantile(Av.temp, 0.25, na.rm = TRUE),
    .groups = "drop"
  )

seasonal_summary$day_index <- 1:nrow(seasonal_summary)

ggplot(seasonal_summary, aes(x = day_index, y = median_temp)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
  xlab("Day of Year") +
  ylab("Median daily temperature (°C)") +
  theme_minimal()

DOCtemp <- DOCtemp %>%
  left_join(seasonal_summary, by = "month_day") %>%
  mutate(
    deseasoned = Av.temp - median_temp
  )

ggplot(DOCtemp, aes(x = Date)) +
  geom_line(aes(y = Av.temp, colour = "Original"), linewidth = 0.3) +
  geom_line(aes(y = deseasoned, colour = "Deseasoned"), linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(colour = "Series", title = "Deseasoned Time Series") +
  theme_minimal()

detrendseasoned <- data.frame(
  original_var = var(DOCtemp$Av.temp, na.rm = TRUE),
  detrended_var = var(DOCtemp$linear_residual, na.rm = TRUE),
  deseasoned_var = var(DOCtemp$deseasoned, na.rm = TRUE)
)

detrendseasoned

# Perform Fourier Transform
fft_result <- fft(DOCtemp$Av.temp)

# Calculate the power spectrum (magnitude squared of Fourier coefficients)
power_spectrum <- data.frame(
  power = Mod(fft_result)^2,
  freq = 1:length(fft_result)
)

ggplot(power_spectrum[1:1000, ], aes(x = freq, y = power)) +
  geom_line() +
  xlab("Fourier component") +
  ylab("Power") +
  labs(title = "Fourier Power Spectrum")+
  theme_minimal()

#auto-correlation and partical autocorrleation functions
acf(DOCtemp$deseasoned, lag.max = 1000)
pacf(DOCtemp$deseasoned, lag.max = 1000)

# Model 2: ARIMA on deseasoned series
deseasoned_ts <- ts(DOCtemp$deseasoned, frequency = 365)
model_arima_deseasoned <- auto.arima(deseasoned_ts)

summary(model_arima_deseasoned)

checkresiduals(model_arima_deseasoned)

DOCtemp$arima_resid <- residuals(model_arima_deseasoned)

ggplot(DOCtemp, aes(x = Date)) +
  geom_line(aes(y = Av.temp, colour = "Original"), linewidth = 0.3) +
  geom_line(aes(y = arima_resid, colour = "ARIMA residual"), linewidth = 0.3) +
  xlab("Date") +
  ylab("Daily average temperature (°C)") +
  scale_colour_manual(values = c("Original" = "firebrick", "ARIMA residual" = "darkblue")) +
  labs(colour = "Series") +
  theme_minimal()

# Forecast deseasoned component for 2020
forecast_arima_2020 <- forecast(model_arima_deseasoned, h = 366)

# Add back seasonal cycle for 2020
seasonal_2020 <- seasonal_summary %>%
  arrange(day_index) %>%
  select(month_day, median_temp)

# Model 3: Fourier + ARIMA

# model selection for K
fourier_models <- lapply(1:8, function(k) {
 auto.arima(temp_ts, xreg = fourier(temp_ts, K = k), seasonal = FALSE)
})

aic_vals <- sapply(fourier_models, AIC)
aic_vals

best_k <- which.min(aic_vals)
best_k   # Previously selected using AIC: best_k = 2

#best_k <- 2

model_fourier <- auto.arima(
  temp_ts,
  xreg = fourier(temp_ts, K = best_k),
  seasonal = FALSE
)

future_fourier <- fourier(temp_ts, K = best_k, h = 366)
forecast_fourier_2020 <- forecast(model_fourier, xreg = future_fourier)

#Load observed data 2020
# https://durhamweather.webspace.durham.ac.uk/2020-data/

temp2020 <- read.csv("durhamtemp_2020.csv", sep=";")

# Convert date
temp2020$Date <- dmy(temp2020$Date)
sum(is.na(temp2020$Date))
head(temp2020$Date)

# Sort
temp2020 <- temp2020[order(temp2020$Date), ]

# Missing dates
full_dates <- seq.Date(min(temp2020$Date), max(temp2020$Date), by = "day")
missing_dates <- setdiff(full_dates, temp2020$Date)
n_missing_dates <- length(missing_dates)

n_missing_dates
head(missing_dates)

# Duplicate dates
n_duplicate_dates <- sum(duplicated(temp2020$Date))
n_duplicate_dates

# Missing values
na_summary <- colSums(is.na(temp2020))
na_summary

# Fill missing values with average of previous and next day
i_max <- which(is.na(temp2020$Max.Temp))
temp2020$Max.Temp[i_max] <- mean(c(temp2020$Max.Temp[i_max - 1],
                               temp2020$Max.Temp[i_max + 1]))

i_min <- which(is.na(temp2020$Min.Temp))
temp2020$Min.Temp[i_min] <- mean(c(temp2020$Min.Temp[i_min - 1],
                               temp2020$Min.Temp[i_min + 1]))
# Recalculate average temperature
temp2020$Av.temp <- (temp2020$Max.Temp + temp2020$Min.Temp) / 2


# Missing values
na_summary <- colSums(is.na(temp2020))
na_summary

nrow(temp2020)

# Create month-day key for seasonal reconstruction
temp2020$month_day <- format(temp2020$Date, "%m-%d")

# 2020 predictions

# Model 1: Seasonal naive
pred_snaive <- as.numeric(model_snaive$mean)

# Model 2: ARIMA on deseasoned series + add seasonal component back
temp2020 <- temp2020 %>%
  left_join(seasonal_summary %>% select(month_day, median_temp), by = "month_day")

pred_arima_deseasoned <- as.numeric(forecast_arima_2020$mean) + temp2020$median_temp

# Model 3: Fourier + ARIMA
pred_fourier <- as.numeric(forecast_fourier_2020$mean)

# Check lengths
length(pred_snaive)
length(pred_arima_deseasoned)
length(pred_fourier)

# Attach predictions
temp2020$Pred_SNaive <- pred_snaive
temp2020$Pred_ARIMA_Deseasoned <- pred_arima_deseasoned
temp2020$Pred_Fourier <- pred_fourier

#Plots observed vs predicted

# Seasonal naive
c1 <- ggplot(temp2020, aes(x = Date)) +
  geom_line(aes(y = Av.temp, colour = "Observed"), linewidth = 0.6) +
  geom_line(aes(y = Pred_SNaive, colour = "Seasonal naive"), linewidth = 0.6) +
  xlab("Date") +
  ylab("Daily average temperature (°C)") +
  scale_colour_manual(values = c("Observed" = "black", "Seasonal naive" = "brown")) +
  labs(colour = "Series", title = "(C) Seasonal Naive Series") +
  theme_minimal()

# ARIMA deseasoned
c2 <- ggplot(temp2020, aes(x = Date)) +
  geom_line(aes(y = Av.temp, colour = "Observed"), linewidth = 0.6) +
  geom_line(aes(y = Pred_ARIMA_Deseasoned, colour = "ARIMA deseasoned"), linewidth = 0.6) +
  xlab("Date") +
  ylab("Daily average temperature (°C)") +
  scale_colour_manual(values = c("Observed" = "black", "ARIMA deseasoned" = "darkgreen")) +
  labs(colour = "Series", title = "(A) ARIMA Deseasoned Series") +
  theme_minimal()

# Fourier + ARIMA
c3 <- ggplot(temp2020, aes(x = Date)) +
  geom_line(aes(y = Av.temp, colour = "Observed"), linewidth = 0.6) +
  geom_line(aes(y = Pred_Fourier, colour = "Fourier + ARIMA"), linewidth = 0.6) +
  xlab("Date") +
  ylab("Daily average temperature (°C)") +
  scale_colour_manual(values = c("Observed" = "black", "Fourier + ARIMA" = "purple")) +
  labs(colour = "Series", title = "(B) Fourier + ARIMA Series") +
  theme_minimal()

grid.arrange(c2, c3, c1, nrow = 3)

# Model comparisons

lm_model_ARIMA_Deseasoned <- lm(temp2020$Pred_ARIMA_Deseasoned ~ temp2020$Av.temp)
lm_model_Fourier_ARIMA <- lm(temp2020$Pred_Fourier ~ temp2020$Av.temp)
lm_model_Fourier_SNaive <- lm(temp2020$Pred_SNaive ~ temp2020$Av.temp)

model_comparison <- data.frame(
  Model = c("Seasonal_Naive", "ARIMA_deseasoned", "Fourier_ARIMA"),
  RMSE = c(
    rmse(temp2020$Av.temp, temp2020$Pred_SNaive),
    rmse(temp2020$Av.temp, temp2020$Pred_ARIMA_Deseasoned),
    rmse(temp2020$Av.temp, temp2020$Pred_Fourier)
  ),
  MAE = c(
    mae(temp2020$Av.temp, temp2020$Pred_SNaive),
    mae(temp2020$Av.temp, temp2020$Pred_ARIMA_Deseasoned),
    mae(temp2020$Av.temp, temp2020$Pred_Fourier)
  ),
  R2 = c(
    summary(lm_model_Fourier_SNaive)$r.squared,
    summary(lm_model_ARIMA_Deseasoned)$r.squared,
    summary(lm_model_Fourier_ARIMA)$r.squared
  )
)

model_comparison <- model_comparison %>%
  arrange(RMSE)

model_comparison

# Best model scatter plot and linear fit
best_model_name <- model_comparison$Model[1]
best_model_name

if (best_model_name == "Seasonal_Naive") {
  best_pred <- temp2020$Pred_SNaive
} else if (best_model_name == "ARIMA_deseasoned") {
  best_pred <- temp2020$Pred_ARIMA_Deseasoned
} else {
  best_pred <- temp2020$Pred_Fourier
}

ggplot(data.frame(Observed = temp2020$Av.temp, Predicted = best_pred),
       aes(x = Observed, y = Predicted)) +
  geom_point(shape = 21, fill = "firebrick") +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, colour = "black") +
  xlab("Observed") +
  ylab("Predicted") +
  labs(title = "ARIMA Deseasoned Model: Scatter-plot of Observed vs Predicted Values") +
  theme_minimal()

model_eval <- lm(best_pred ~ temp2020$Av.temp)
summary(model_eval)
