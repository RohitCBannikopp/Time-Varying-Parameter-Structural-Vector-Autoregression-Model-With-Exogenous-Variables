#############################
# TVP-SVARX Post-processing
# Requires: `chains1` (sampler output) already in memory
# Outputs:
#   - Time-varying FEVDs (domestic vs global) at all horizons 0..H
#   - Period summaries + P(global > domestic)
#   - IRF summaries by event windows
#   - Plots (FEVD time series, stacked area, period comparison, IRF ribbons)
#############################

# =============== Libraries ===============
library(abind)
library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(progress)


# =============== SETTINGS ===============
L <- 3
var_names <- c("N_50","real_rate","term_spread","Dp_ratio","rel_bill_rate")
n <- length(var_names)
target_var <- "N_50"
target_idx <- which(var_names == target_var)

exo_names <- c("GDP","IIP","CPI","D_monetary","US_monetary","Oil_s","Oil_d")
m_exo <- length(exo_names)

# Alphahat index mapping (per-equation state vector)
intercept_idx <- 1
lag_start <- 2
lag_end   <- lag_start + n*L - 1     # 2:16 for n=5, L=3
exo_start <- lag_end + 1             # 17
exo_end   <- exo_start + m_exo - 1   # 17:23 for m_exo=7

# Sample dates (edit if your sample differs)
dates <- seq(as.Date("2005-11-01"), as.Date("2024-03-01"), by = "month")
Tt <- length(dates)
if (Tt != 221) warning("Check your `dates`; expected 221 months.")

# Posterior settings (must match your sampler)
burnin <- 3000
thin_draws <- 10
ci_lower <- 0.05
ci_upper <- 0.95

# --------- Horizon conventions (IMPORTANT) ----------
# H = simulation horizon used inside posterior calculations for IRFs/FEVDs.
#    We simulate responses 0..H each iteration. This sets how far ahead
#    the model’s dynamic responses are computed per draw. (Used in loops.)
H <- 24
# choose_h = interpretation horizon. After we have FEVD/IRF arrays for all 0..H,
#    we summarize and plot at this single horizon (e.g., 6 months ahead).
choose_h <- 6
h_idx <- choose_h + 1
# ---------------------------------------------------

# Mapping exogenous shocks to Domestic vs Global
domestic_exo_cols <- 1:4             # GDP, IIP, CPI, D_monetary
global_exo_cols   <- 5:7             # US_monetary, Oil_s, Oil_d

# Event windows
event_windows <- list(
  pre_GFC        = c(as.Date("2005-11-01"), as.Date("2007-07-31")),
  GFC            = c(as.Date("2007-08-01"), as.Date("2009-03-31")),
  post_GFC       = c(as.Date("2009-04-01"), as.Date("2013-05-31")),
  taper_tantrum  = c(as.Date("2013-06-01"), as.Date("2013-12-31")),
  pre_COVID      = c(as.Date("2014-01-01"), as.Date("2019-12-31")),
  COVID          = c(as.Date("2020-01-01"), as.Date("2020-12-31")),
  post_COVID     = c(as.Date("2021-01-01"), as.Date("2024-03-01"))
)
shocks_of_interest_exo <- 1:m_exo

# =============== Helpers ===============
build_companion <- function(Phi_block, n, L){
  k <- n * L
  Acomp <- matrix(0, n*L, n*L)
  Acomp[1:n, 1:(n*L)] <- Phi_block
  if (L > 1) Acomp[(n+1):(n*L), 1:(n*(L-1))] <- diag(n*(L-1))
  Acomp
}

compute_irf_endog <- function(Acomp, S, H){
  n_local <- nrow(S)
  k <- nrow(Acomp)
  J <- cbind(diag(n_local), matrix(0, n_local, k - n_local))
  S_big <- rbind(S, matrix(0, k - n_local, n_local))
  irf <- array(0, c(n_local, n_local, H+1))
  irf[,,1] <- S
  A_power <- diag(k)
  for (h in 1:H){
    A_power <- A_power %*% Acomp
    irf[,,h+1] <- J %*% A_power %*% S_big
  }
  irf
}

compute_irf_exog <- function(Acomp, B, shock_col_index, H){
  n_local <- nrow(B)
  k <- nrow(Acomp)
  J <- cbind(diag(n_local), matrix(0, n_local, k - n_local))
  b_col <- B[, shock_col_index, drop = FALSE]
  b_big <- rbind(b_col, matrix(0, k - n_local, 1))
  irf_exo <- matrix(0, nrow = n_local, ncol = H+1)
  irf_exo[,1] <- J %*% b_big
  A_power <- diag(k)
  for (h in 1:H){
    A_power <- A_power %*% Acomp
    irf_exo[,h+1] <- J %*% A_power %*% b_big
  }
  irf_exo
}

# =============== MAIN POSTERIOR LOOP ===============
agg_global_list  <- list()
agg_domestic_list <- list()
agg_other_list   <- list()

# store event-averaged IRFs: list[event] -> array(draws x (H+1) x shocks)
event_irf_draws <- setNames(vector("list", length(event_windows)), names(event_windows))

draw_counter <- 0
total_draws <- length(chains1) * (length(chains1[[1]]$theta) - burnin) / thin_draws


# Setup progress bar across all posterior draws
total_iters <- sum(sapply(chains1, function(ch) {
  length(seq(burnin + 1, length(ch$theta), by = thin_draws))
}))
pb <- progress_bar$new(
  format = "  Processing [:bar] :percent eta: :eta",
  total = total_draws, clear = FALSE, width = 60
)

for (ch in seq_along(chains1)) {
  ch_main <- chains1[[ch]]
  niter_chain <- length(ch_main$theta)
  iter_seq <- seq(burnin + 1, niter_chain, by = thin_draws)
  
  for (iter in iter_seq) {
    draw_counter <- draw_counter + 1
    
    theta_list <- ch_main$theta[[iter]]   # list length n; each has T x m_state `alphahat`
    u_mat      <- ch_main$u[[iter]]       # T x n residuals
    A_list     <- ch_main$A[[iter]]       # list length T: structural A_t (if available)
    
    # per-draw FEVD storage (H+1 x T)
    fevd_glob <- matrix(NA_real_, H+1, Tt)
    fevd_dom  <- matrix(NA_real_, H+1, Tt)
    fevd_oth  <- matrix(NA_real_, H+1, Tt)
    
    # accumulators for event windows
    event_accum  <- lapply(names(event_windows), function(x) array(0, c(length(shocks_of_interest_exo), H+1)))
    names(event_accum) <- names(event_windows)
    event_counts <- setNames(as.list(rep(0, length(event_windows))), names(event_windows))
    
    # ---- time loop ----
    for (t in 1:Tt) {
      # Build Phi (n x nL) and B (n x m_exo) from alphahat
      Phi_block <- matrix(NA_real_, n, n*L)
      B_mat     <- matrix(0, n, m_exo)
      for (i in 1:n) {
        coeffs_t <- as.numeric(theta_list[[i]]$alphahat[t, ])
        Phi_block[i, ] <- coeffs_t[lag_start:lag_end]
        B_mat[i, ]     <- coeffs_t[exo_start:exo_end]
      }
      Acomp <- build_companion(Phi_block, n, L)
      
      # Structural impact matrix S = A_t^{-1}
      S <- NULL
      if (!is.null(A_list) && length(A_list) >= t && is.matrix(A_list[[t]])) {
        S <- tryCatch(solve(A_list[[t]]), error = function(e) diag(1e-6, n))
      } else {
        S <- diag(1e-6, n)
      }
      
      # IRFs to endogenous & exogenous shocks
      irf_endog   <- compute_irf_endog(Acomp, S, H)              # n x n x (H+1)
      irf_exo_all <- array(0, c(n, H+1, m_exo))
      for (j in 1:m_exo) {
        irf_exo_all[,,j] <- compute_irf_exog(Acomp, B_mat, j, H)
      }
      
      # FEVD for target variable
      tot_shocks <- n + m_exo
      cum_sq <- matrix(0, tot_shocks, H+1)
      
      # endogenous shocks
      for (sh in 1:n) {
        for (h in 0:H) {
          cum_sq[sh, h+1] <- sum(irf_endog[target_idx, sh, 1:(h+1)]^2)
        }
      }
      # exogenous shocks
      for (j in 1:m_exo) {
        for (h in 0:H) {
          cum_sq[n + j, h+1] <- sum(irf_exo_all[target_idx, 1:(h+1), j]^2)
        }
      }
      
      total_var <- colSums(cum_sq)
      shares <- sweep(cum_sq, 2, total_var, "/")   # shock-by-horizon shares
      
      dom_share  <- colSums(shares[n + domestic_exo_cols, , drop = FALSE], na.rm = TRUE)
      glob_share <- colSums(shares[n + global_exo_cols,   , drop = FALSE], na.rm = TRUE)
      oth_share  <- pmax(0, 1 - glob_share - dom_share)
      
      fevd_dom[,t]  <- dom_share
      fevd_glob[,t] <- glob_share
      fevd_oth[,t]  <- oth_share
      
      # Accumulate event-averaged IRFs (for target var)
      date_t <- dates[t]
      for (ev in names(event_windows)) {
        rng <- event_windows[[ev]]
        if (date_t >= rng[1] && date_t <= rng[2]) {
          event_counts[[ev]] <- event_counts[[ev]] + 1
          for (si in seq_along(shocks_of_interest_exo)) {
            j <- shocks_of_interest_exo[si]
            event_accum[[ev]][si, ] <- event_accum[[ev]][si, ] + irf_exo_all[target_idx, , j]
          }
        }
      }
    } # end time loop
    
    # Save FEVD results
    agg_global_list[[length(agg_global_list) + 1]]  <- fevd_glob
    agg_domestic_list[[length(agg_domestic_list) + 1]] <- fevd_dom
    agg_other_list[[length(agg_other_list) + 1]]   <- fevd_oth
    
    # Save event-averaged IRFs
    for (ev in names(event_windows)) {
      cnt <- event_counts[[ev]]
      if (cnt > 0) {
        avg_irf_mat <- event_accum[[ev]] / cnt   # shocks x (H+1)
        new_draw_array <- array(NA_real_, c(1, H+1, nrow(avg_irf_mat)))
        for (si in 1:nrow(avg_irf_mat)) new_draw_array[1, , si] <- avg_irf_mat[si, ]
        if (is.null(event_irf_draws[[ev]])) {
          event_irf_draws[[ev]] <- new_draw_array
        } else {
          event_irf_draws[[ev]] <- abind(event_irf_draws[[ev]], new_draw_array, along = 1)
        }
      }
    }
    

    # only update every 50 draws
    if (draw_counter == 1 || draw_counter %% 50 == 0) {
      pb$update(draw_counter / total_draws)
    }
    
  }
}



ndraws_kept <- length(agg_global_list)
message("Posterior draws saved: ", ndraws_kept)

# =============== Summaries at interpretation horizon h = choose_h = ", choose_h, " ===============
# Convert lists (each (H+1 x T)) into arrays: (H+1 x T x draws)
agg_global_array  <- simplify2array(agg_global_list)
agg_domestic_array <- simplify2array(agg_domestic_list)
# agg_other_array  <- simplify2array(agg_other_list)  # available if you need it

# Time series at h = choose_h
global_median_ts <- apply(agg_global_array[h_idx, , ], 1, median, na.rm = TRUE)
global_low_ts    <- apply(agg_global_array[h_idx, , ], 1, quantile, probs = ci_lower, na.rm = TRUE)
global_high_ts   <- apply(agg_global_array[h_idx, , ], 1, quantile, probs = ci_upper, na.rm = TRUE)

dom_median_ts <- apply(agg_domestic_array[h_idx, , ], 1, median, na.rm = TRUE)
dom_low_ts    <- apply(agg_domestic_array[h_idx, , ], 1, quantile, probs = ci_lower, na.rm = TRUE)
dom_high_ts   <- apply(agg_domestic_array[h_idx, , ], 1, quantile, probs = ci_upper, na.rm = TRUE)

fevd_df <- data.frame(
  date = dates,
  global_median = global_median_ts, global_low = global_low_ts, global_high = global_high_ts,
  dom_median = dom_median_ts,       dom_low = dom_low_ts,       dom_high = dom_high_ts
)

# Period summaries (means across dates, distribution across draws)
periods <- event_windows
period_summary <- data.frame(
  period = names(periods),
  start  = as.Date(sapply(periods, `[`, 1)),
  end    = as.Date(sapply(periods, `[`, 2)),
  global_mean = NA_real_, global_ci_low = NA_real_, global_ci_high = NA_real_,
  domestic_mean = NA_real_, dom_ci_low = NA_real_, dom_ci_high = NA_real_,
  P_global_gt_domestic = NA_real_
)

for (i in seq_along(periods)){
  d0 <- periods[[i]][1]; d1 <- periods[[i]][2]
  idx <- which(dates >= d0 & dates <= d1)
  if (length(idx) == 0) next
  
  # average across time points in idx but retain cross-draw distribution
  glob_draws_period <- apply(agg_global_array[h_idx, idx, ], 2, mean, na.rm = TRUE)
  dom_draws_period  <- apply(agg_domestic_array[h_idx, idx, ], 2, mean, na.rm = TRUE)
  
  period_summary$global_mean[i]    <- median(glob_draws_period, na.rm = TRUE)
  period_summary$global_ci_low[i]  <- quantile(glob_draws_period, probs = ci_lower, na.rm = TRUE)
  period_summary$global_ci_high[i] <- quantile(glob_draws_period, probs = ci_upper, na.rm = TRUE)
  
  period_summary$domestic_mean[i]  <- median(dom_draws_period, na.rm = TRUE)
  period_summary$dom_ci_low[i]     <- quantile(dom_draws_period, probs = ci_lower, na.rm = TRUE)
  period_summary$dom_ci_high[i]    <- quantile(dom_draws_period, probs = ci_upper, na.rm = TRUE)
  
  period_summary$P_global_gt_domestic[i] <- mean(glob_draws_period > dom_draws_period, na.rm = TRUE)
}

print(period_summary)
write_xlsx(period_summary, path = "FEVD Data.xlsx")

# =============== Plots ===============
# FEVD time series at h = choose_h
fevd_long <- fevd_df |>
  select(date,
         global_median, global_low, global_high,
         dom_median, dom_low, dom_high) |>
  pivot_longer(cols = -date, names_to = c("type","stat"), names_sep = "_", values_to = "value") |>
  pivot_wider(names_from = stat, values_from = value)

ggplot(fevd_long, aes(x = date, color = type, fill = type)) +
  geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.15, linetype = 0) +
  geom_line(aes(y = median), linewidth = 0.8) +
  scale_color_manual(values = c("global" = "#3B82F6", "dom" = "#F97316")) +
  scale_fill_manual(values  = c("global" = "#3B82F6", "dom" = "#F97316")) +
  labs(title = paste0("FEVD shares for ", target_var, " (", choose_h, "-month horizon)"),
       y = "Share of forecast variance", x = "Date",
       color = "Shock Type", fill = "Shock Type") +
  theme_minimal()

# Stacked area (median)
hist_decomp_df <- data.frame(
  date = rep(dates, 3),
  share = c(global_median_ts, dom_median_ts, pmax(0, 1 - global_median_ts - dom_median_ts)),
  component = rep(c("Global","Domestic","Other"), each = length(dates))
)
ggplot(hist_decomp_df, aes(x = date, y = share, fill = component)) +
  geom_area() +
  theme_minimal() +
  labs(title = paste("Median FEVD decomposition for", target_var, "(h =", choose_h, ")"),
       y = "Share", x = "Date")

# Period plot with CIs (uses `period_summary`)
plot_data <- period_summary |>
  transmute(period,
            type = "Global Shocks", mean = global_mean, ci_low = global_ci_low, ci_high = global_ci_high) |>
  bind_rows(
    period_summary |>
      transmute(period,
                type = "Domestic Shocks", mean = domestic_mean, ci_low = dom_ci_low, ci_high = dom_ci_high)
  )

ggplot(plot_data, aes(x = period, y = mean, color = type, group = type)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.2, position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Global Shocks" = "steelblue", "Domestic Shocks" = "darkred")) +
  theme_minimal(base_size = 14) +
  labs(title = "Impact of Global vs Domestic Shocks Across Periods",
       x = "Period", y = "Posterior Mean with 90% Credible Interval",
       color = "Shock Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(size = 16, face = "bold", hjust = 0.5))

# =============== IRF ANALYSIS (event-averaged, exogenous shocks) ===============
# IRF summary table at horizon h_target (e.g., 6)
h_target <- 6
irf6_summary_list <- list()

for (ev in names(event_irf_draws)){
  mat_ev <- event_irf_draws[[ev]]   # draws x (H+1) x shocks
  if (is.null(mat_ev)) next
  n_shocks_ev <- dim(mat_ev)[3]
  
  ev_rows <- lapply(seq_len(n_shocks_ev), function(sind){
    draws_h_mat <- mat_ev[ , , sind]              # draws x (H+1)
    h_draws <- draws_h_mat[, h_target + 1]        # IRF at horizon = h_target
    data.frame(
      event       = ev,
      shock       = exo_names[shocks_of_interest_exo[sind]],
      irf_median  = median(h_draws, na.rm = TRUE),
      irf_low     = quantile(h_draws, probs = ci_lower, na.rm = TRUE),
      irf_high    = quantile(h_draws, probs = ci_upper, na.rm = TRUE),
      row.names = NULL, check.names = FALSE
    )
  })
  
  irf6_summary_list[[ev]] <- bind_rows(ev_rows)
}

irf6_table <- bind_rows(irf6_summary_list)
print(irf6_table)
write_xlsx(irf6_table, path = "IRF data.xlsx")

# IRF ribbon plots (0..max_h_plot) by event
max_h_plot <- 6
for (ev in names(event_irf_draws)) {
  mat_ev <- event_irf_draws[[ev]]
  if (is.null(mat_ev)) next
  
  plot_df <- bind_rows(lapply(seq_len(dim(mat_ev)[3]), function(sind){
    dh <- mat_ev[,,sind]  # draws x (H+1)
    data.frame(
      h      = 0:max_h_plot,
      median = apply(dh, 2, median, na.rm = TRUE)[1:(max_h_plot + 1)],
      low    = apply(dh, 2, quantile, probs = ci_lower, na.rm = TRUE)[1:(max_h_plot + 1)],
      high   = apply(dh, 2, quantile, probs = ci_upper, na.rm = TRUE)[1:(max_h_plot + 1)],
      shock  = exo_names[shocks_of_interest_exo[sind]]
    )
  }))
  
  p <- ggplot(plot_df, aes(x = h, y = median, color = shock, fill = shock)) +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = 0:max_h_plot) +
    labs(title = paste("IRFs –", ev, "(0–", max_h_plot, " Month Horizon)"),
         x = "Horizon (months)", y = paste("Response of", target_var),
         color = "Shock", fill = "Shock") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  print(p)
}

# =============== End ===============
