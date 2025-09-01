# Load necessary libraries
install.packages("tseries")  # For adf.test, pp.test
install.packages("urca")     # For ur.kpss
library(tseries)
library(urca)


# Test for each matrix (endo_clean and exo_clean)
library(tseries)

apply(exo_clean, 2, function(x) {
  adf.test(x, alternative = "stationary")
})

apply(exo_clean, 2, function(x) {
  kpss.test(x, null = "Level")
})
apply(endo_clean, 2, function(x) {
  adf.test(x, alternative = "stationary")
})

apply(endo_clean, 2, function(x) {
  kpss.test(x, null = "Level")
})
