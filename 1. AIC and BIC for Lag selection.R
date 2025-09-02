# ---------------------------------------------
# Lag Length Selection using Information Criteria
# ---------------------------------------------
library(vars)

# 'endo_clean' = T x n matrix or dataframe of endogenous variables (after cleaning lags)

# Select optimal lag length up to 6 (you can increase if needed)
lag_selection <- VARselect(endo_clean, lag.max = 6, type = "const")

# Print full table of results (shows AIC, HQ, SC, FPE for each lag order)
print(lag_selection)

# Extract chosen lag length according to different criteria
lag_AIC <- lag_selection$selection["AIC(n)"]   # Akaike Information Criterion
lag_BIC <- lag_selection$selection["SC(n)"]    # Schwarz Bayesian Information Criterion (BIC)

# Display results
cat("Optimal lag suggested by AIC:", lag_AIC, "\n")
cat("Optimal lag suggested by BIC:", lag_BIC, "\n")

# Notes:
# - AIC often selects a higher lag order (better fit, less parsimonious).
# - BIC tends to select a smaller lag order (more parsimonious).
# - Final choice depends on balancing fit vs. degrees of freedom.
# ---------------------------------------------
