# =======================================================
# Posterior Residual Diagnostics (ACF + Ljung–Box test)
# =======================================================
library(dplyr)
library(ggplot2)

# ----- Settings -----
max_lag   <- 6          # check autocorrelations up to lag 6
ci_lower  <- 0.05       # lower quantile for posterior credible band
ci_upper  <- 0.95       # upper quantile for posterior credible band

# Function: compute sample autocorrelation at lag k
acf_at_lag <- function(x, k){
  x <- as.numeric(x) - mean(x, na.rm = TRUE)
  n <- length(x)
  if(k >= n) return(NA_real_)
  num <- sum(x[(k+1):n] * x[1:(n-k)])
  den <- sum(x^2)
  if(den == 0) return(NA_real_)
  num / den
}

# ----- Loop over all endogenous variables -----
acf_summary_list <- list()
acf_plot_list    <- list()

for (v in seq_along(var_names)) {
  varname <- var_names[v]
  
  acf_draws_mat <- NULL   # store posterior draws of ACFs
  draw_ids <- c()
  
  # Loop across chains
  for (ch in seq_along(chains1)) {
    ch_raw  <- chains1[[ch]]
    ch_main <- if (!is.null(ch_raw$theta)) ch_raw else ch_raw[[1]]
    niter_chain <- length(ch_main$u)
    
    # Only keep iterations after burnin and apply thinning
    iter_seq <- seq(burnin + 1, niter_chain, by = thin_draws)
    
    for (iter in iter_seq) {
      u_mat <- ch_main$u[[iter]]   # residuals (T x n)
      
      # Ensure correct orientation (rows = time, cols = variables)
      if (ncol(u_mat) != length(var_names)) {
        if (nrow(u_mat) == length(var_names) && ncol(u_mat) == Tt) {
          u_mat <- t(u_mat)
        }
      }
      if (nrow(u_mat) != Tt) next
      
      # Extract residuals for this variable
      u_t <- u_mat[, v]
      
      # Compute autocorrelations up to max_lag
      vals <- vapply(1:max_lag, function(L) acf_at_lag(u_t, L), numeric(1))
      
      # Collect across posterior draws
      acf_draws_mat <- rbind(acf_draws_mat, vals)
      draw_ids <- c(draw_ids, paste0("c", ch, "_i", iter))
    }
  }
  
  n_draws <- nrow(acf_draws_mat)
  message("Variable ", varname, " | Posterior draws kept for ACF: ", n_draws)
  
  # Summarize posterior ACF distribution
  lags <- 1:max_lag
  acf_summary <- data.frame(
    lag      = lags,
    median   = apply(acf_draws_mat, 2, median, na.rm = TRUE),
    low      = apply(acf_draws_mat, 2, quantile, probs = ci_lower, na.rm = TRUE),
    high     = apply(acf_draws_mat, 2, quantile, probs = ci_upper, na.rm = TRUE),
    variable = varname
  )
  
  acf_summary_list[[varname]] <- acf_summary
  
  # ----- Plot ACF with posterior bands -----
  se_band <- 1.96 / sqrt(Tt)   # white-noise CI band
  p_acf <- ggplot(acf_summary, aes(x = lag, y = median)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.18, fill = "skyblue") +
    geom_line(linewidth = 0.9, color = "blue") +
    geom_point(size = 2, color = "blue") +
    geom_hline(yintercept =  se_band, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -se_band, linetype = "dashed", linewidth = 0.3) +
    labs(title    = paste("Posterior ACF of residuals (", varname, ")", sep=""),
         subtitle = "Median and 90% credible bands across TVP-SVARX posterior draws",
         x = "Lag", y = "Autocorrelation") +
    theme_minimal()
  
  acf_plot_list[[varname]] <- p_acf
  
  # ----- Ljung–Box Test (posterior distribution) -----
  Q_vals <- apply(acf_draws_mat, 1, function(rhos){
    k <- 1:max_lag
    Tt * (Tt + 2) * sum((rhos^2) / (Tt - k))
  })
  
  p_vals <- 1 - pchisq(Q_vals, df = max_lag)
  post_prob_reject_5pct <- mean(p_vals < 0.05, na.rm = TRUE)
  
  cat(sprintf(
    "\n%s | Posterior Ljung–Box(%d): median Q = %.3f, median p = %.3f\nPosterior Pr(p<0.05) = %.2f\n",
    varname, max_lag, median(Q_vals, na.rm = TRUE), median(p_vals, na.rm = TRUE), post_prob_reject_5pct
  ))
}

# ----- Example: show plots for all variables -----
for (v in var_names) {
  print(acf_plot_list[[v]])
}
