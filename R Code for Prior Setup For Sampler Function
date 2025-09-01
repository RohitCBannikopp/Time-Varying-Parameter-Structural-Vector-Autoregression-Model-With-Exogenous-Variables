library(KFAS)

# Dimensions
T <- nrow(endo_clean)
n <- ncol(endo_clean)
m <- ncol(exo_clean)

# Regressor matrix: intercept + lags + exogenous variables
intercept <- matrix(1, nrow = T, ncol = 1)
Z_mat <- cbind(intercept, lag_endo_clean, exo_clean)
q <- ncol(Z_mat)

# State priors
a1 <- rep(0, q)              # prior mean of coefficients
P1 <- diag(1e6, q)           # diffuse prior covariance
Q  <- diag(1e-4, q)          # innovation variance (fixed across equations)

# Build state-space models for each endogenous variable
models <- vector("list", n)
for (i in 1:n) {
  y_i <- endo_clean[, i]
  
  Z_i <- array(NA, dim = c(1, q, T))
  for (t in 1:T) {
    Z_i[1, , t] <- Z_mat[t, ]
  }
  
  models[[i]] <- SSModel(
    y_i ~ -1 + SSMcustom(
      Z = Z_i,
      T = diag(q),
      Q = Q,
      a1 = a1,
      P1 = P1
    ),
    H = 1
  )
}
