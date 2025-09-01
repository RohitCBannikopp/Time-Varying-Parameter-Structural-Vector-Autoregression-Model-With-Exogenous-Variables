library(readxl)
library(dplyr)

# Load data
Vdata <- read_excel("C:/Users/HP/Desktop/TVP-SVARX modelling/Variables.xlsx")

# Endogenous and exogenous variables
endo <- as.matrix(Vdata[, c("N_50","real_rate","term_spread","Dp_ratio","re_bill_rate")])
exo  <- as.matrix(Vdata[, c("GDP","IIP","CPI","D_monetory","US_monetory","Oil_s","Oil_d")])

# Number of lags
k <- 3  

# Function to create lagged matrices
create_lagged_matrix <- function(data, lags) {
  lagged_list <- lapply(1:lags, function(lag_i) {
    dplyr::lag(data, n = lag_i) %>%
      setNames(paste0(colnames(data), "_lag", lag_i))
  })
  bind_cols(lagged_list)
}

# Create lagged matrices
lag_endo <- create_lagged_matrix(as.data.frame(endo), k)
lag_exo  <- create_lagged_matrix(as.data.frame(exo), k)

# Drop initial rows lost to lagging
lag_endo_clean <- lag_endo[-(1:k), , drop = FALSE]
lag_exo_clean  <- lag_exo[-(1:k), , drop = FALSE]
endo_clean     <- endo[-(1:k), , drop = FALSE]
exo_clean      <- exo[-(1:k), , drop = FALSE]
