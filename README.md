# Daily Temperature Forecasting in Durham Using Classical Time Series Models

## Overview

This project forecasts daily mean temperature in Durham using historical observations from 1901–2019 and evaluates predictions against observed data from 2020.

The objective is to assess how well classical time series models capture long-term trends, strong annual seasonality, and short-term variability in environmental data.

## Data

* Source: Durham Observatory
* Period:

  * Training: 1901–2019
  * Testing: 2020
* Frequency: Daily

## Methods

The following models are compared:

* Seasonal Naive
* ARIMA on deseasonalised series
* Fourier terms with ARIMA errors

## Evaluation

Model performance is evaluated using:

* RMSE
* MAE

## Key Results (to be updated)

* Best performing model:
* RMSE:
* MAE:

## Tools

* R
* tidyverse
* forecast
* tseries
* zoo
* ggplot2

## Project Structure

* `data/` → raw and processed datasets
* `scripts/` → analysis workflow
* `outputs/` → figures and tables
* `report/` → final report

## Status

Finished
