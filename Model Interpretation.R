#############################
# TVP-SVARX Post-processing
# Requires: chains1 (sampler output) in memory
# Outputs:
#   - Time-varying FEVDs (domestic vs global)
#   - Period summaries with posterior probabilities
#   - Plots: FEVD time series, stacked decomposition, event summaries
#############################

# ================== Libraries ==================
library(abind)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(writexl)

# ================== SETTINGS ==================
L <- 3
var_names <- c("N_50","real_rate","term_spread","Dp_ratio","rel_bill_rate")
n <- length(var_names)
target_var <- "N_50"
target_idx <- which(var_names == target_var)

exo_names <- c("GDP","IIP","CPI","D_monetary","US_monetary","Oil_s","Oil_d")
m_exo <- length(exo_names)

# Alphahat index mapping
intercept_idx <- 1
lag_start <- 2
lag_end   <- lag_start + n*L - 1     # 2:16
exo_start <- lag_end + 1             # 17
exo_end   <- exo_start + m_exo - 1   # 17:23

# Sample dates
dates <- seq(as.Date("2005-11-01"), as.Date("2024-03-01"), by = "month")
Tt <- length(dates)
if(Tt != 221) warning("Check your dates vector; expected 221 months.")

# Posterior settings
burnin <- 3000
thin_draws <- 10
ci_lower <- 0.05
ci_upper <- 0.95

# FEVD/IRF settings
# H = 24 is used here for the PRIOR / POSTERIOR SAMPLER setup
# (state-space system horizon in estimation).
# ⚠️ Do not confuse this H with 'choose_h' used later in FEVD/IRF interpretation.
H <- 24
# ----------------------------------------------------------

# ... [your estimation code continues here] ...

# Later in FEVD / IRF analysis:
# choose_h specifies the FORECAST HORIZON (h) at which the variance 
# decomposition or impulse responses are interpreted.
# Example:
choose_h <- 6

# Mapping shocks
domestic_exo_cols <- 1:4
global_exo_cols   <- 5:7

event_windows <- list(
  pre_GFC    = c(as.Date("2005-11-01"), as.Date("2007-07-31")),
  GFC        = c(as.Date("2007-08-01"), as.Date("2009-03-31")),
  post_GFC   = c(as.Date("2009-04-01"), as.Date("2013-05-31")),
  taper_tantrum = c(as.Date("2013-06-01"), as.Date("2013-12-31")),
  pre_COVID  = c(as.Date("2014-01-01"), as.Date("2019-12-31")),
  COVID      = c(as.Date("2020-01-01"), as.Date("2020-12-31")),
  post_COVID = c(as.Date("2021-01-01"), as.Date("2024-03-01"))
)

shocks_of_interest_exo <- 1:m_exo

# ================== Helper Functions ==================
build_companion <- function(Phi_block, n, L){
  k <- n * L
  Acomp <- matrix(0, n*L, n*L)
  Acomp[1:n, 1:(n*L)] <- Phi_block
  if(L > 1) Acomp[(n+1):(n*L), 1:(n*(L-1))] <- diag(n*(L-1))
  return(Acomp)
}

compute_irf_endog <- function(Acomp, S, H){
  n_local <- nrow(S)
  k <- nrow(Acomp)
  J <- cbind(diag(n_local), matrix(0, n_local, k - n_local))
  S_big <- rbind(S, matrix(0, k - n_local, n_local))
  irf <- array(0, c(n_local, n_local, H+1))
  irf[,,1] <- S
  A_power <- diag(k)
  for(h in 1:H){
    A_power <- A_power %*% Acomp
    irf[,, h+1] <- J %*% A_power %*% S_big
  }
  return(irf)
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
  for(h in 1:H){
    A_power <- A_power %*% Acomp
    irf_exo[, h+1] <- J %*% A_power %*% b_big
  }
  return(irf_exo)
}

# ================== MAIN SAMPLER LOOP ==================
agg_global_list <- list()
agg_domestic_list <- list()
agg_other_list <- list()
event_irf_draws <- setNames(vector("list", length(event_windows)), names(event_windows))

draw_counter <- 0
for(ch in seq_along(chains1)){
  ch_main <- chains1[[ch]]
  niter_chain <- length(ch_main$theta)
  
  iter_seq <- seq(burnin + 1, niter_chain, by = thin_draws)
  for(iter in iter_seq){
    draw_counter <- draw_counter + 1
    theta_list <- ch_main$theta[[iter]]
    u_mat <- ch_main$u[[iter]]
    A_list <- ch_main$A[[iter]]
    
    fevd_global_draw <- matrix(NA_real_, H+1, Tt)
    fevd_domestic_draw <- matrix(NA_real_, H+1, Tt)
    fevd_other_draw <- matrix(NA_real_, H+1, Tt)
    
    # Event accumulators
    event_accum <- lapply(names(event_windows), function(x) array(0, c(length(shocks_of_interest_exo), H+1)))
    names(event_accum) <- names(event_windows)
    event_counts <- setNames(as.list(rep(0, length(event_windows))), names(event_windows))
    
    # ----- Loop over time -----
    for(t in 1:Tt){
      # Phi & B
      Phi_block <- matrix(NA, n, n*L)
      B_mat <- matrix(0, n, m_exo)
      for(i in 1:n){
        coeffs_t <- as.numeric(theta_list[[i]]$alphahat[t, ])
        Phi_block[i, ] <- coeffs_t[lag_start:lag_end]
        B_mat[i, ] <- coeffs_t[exo_start:exo_end]
      }
      
      Acomp <- build_companion(Phi_block, n, L)
      
      # Structural impact matrix
      S <- tryCatch(solve(A_list[[t]]), error = function(e) diag(1e-6, n))
      
      # IRFs
      irf_endog <- compute_irf_endog(Acomp, S, H)
      irf_exo_all <- array(0, c(n, H+1, m_exo))
      for(j in 1:m_exo) irf_exo_all[,,j] <- compute_irf_exog(Acomp, B_mat, j, H)
      
      # FEVD
      tot_shocks <- n + m_exo
      cum_sq <- matrix(0, tot_shocks, H+1)
      for(sh in 1:n){
        for(h in 0:H){
          cum_sq[sh, h+1] <- sum(irf_endog[target_idx, sh, 1:(h+1)]^2)
        }
      }
      for(j in 1:m_exo){
        for(h in 0:H){
          cum_sq[n+j, h+1] <- sum(irf_exo_all[target_idx, 1:(h+1), j]^2)
        }
      }
      total_var <- colSums(cum_sq)
      shares <- sweep(cum_sq, 2, total_var, "/")
      
      dom_share <- colSums(shares[n + domestic_exo_cols, , drop=FALSE], na.rm=TRUE)
      glob_share <- colSums(shares[n + global_exo_cols, , drop=FALSE], na.rm=TRUE)
      other_share <- pmax(0, 1 - glob_share - dom_share)
      
      fevd_domestic_draw[,t] <- dom_share
      fevd_global_draw[,t] <- glob_share
      fevd_other_draw[,t] <- other_share
      
      # Events accumulation
      date_t <- dates[t]
      for(ev in names(event_windows)){
        if(date_t >= event_windows[[ev]][1] & date_t <= event_windows[[ev]][2]){
          event_counts[[ev]] <- event_counts[[ev]] + 1
          for(si in seq_along(shocks_of_interest_exo)){
            exo_j <- shocks_of_interest_exo[si]
            irf_vec <- irf_exo_all[target_idx,,exo_j]
            event_accum[[ev]][si, ] <- event_accum[[ev]][si, ] + irf_vec
          }
        }
      }
    } # end time loop
    
    # Save per-draw FEVDs
    agg_global_list[[length(agg_global_list)+1]] <- fevd_global_draw
    agg_domestic_list[[length(agg_domestic_list)+1]] <- fevd_domestic_draw
    agg_other_list[[length(agg_other_list)+1]] <- fevd_other_draw
    
    # Event IRFs per draw
    for(ev in names(event_windows)){
      cnt <- event_counts[[ev]]
      if(cnt > 0){
        avg_irf_mat <- event_accum[[ev]] / cnt
        new_draw_array <- array(NA_real_, c(1, H+1, nrow(avg_irf_mat)))
        for(si in 1:nrow(avg_irf_mat)) new_draw_array[1,,si] <- avg_irf_mat[si, ]
        event_irf_draws[[ev]] <- abind(event_irf_draws[[ev]], new_draw_array, along=1)
      }
    }
  }
}

ndraws_kept <- length(agg_global_list)
message("Posterior draws saved: ", ndraws_kept)

# ================== Summaries ==================
agg_global_array <- simplify2array(agg_global_list)
agg_domestic_array <- simplify2array(agg_domestic_list)

h_idx <- choose_h + 1
global_median_ts <- apply(agg_global_array[h_idx,,], 1, median, na.rm=TRUE)
dom_median_ts <- apply(agg_domestic_array[h_idx,,], 1, median, na.rm=TRUE)

fevd_df <- data.frame(date=dates,
                      global=global_median_ts,
                      domestic=dom_median_ts)

# Period summaries
period_summary <- data.frame(period=names(event_windows),
                             start=sapply(event_windows, `[`, 1),
                             end=sapply(event_windows, `[`, 2),
                             global_mean=NA, domestic_mean=NA,
                             P_global_gt_domestic=NA)
for(i in seq_along(event_windows)){
  idx <- which(dates >= event_windows[[i]][1] & dates <= event_windows[[i]][2])
  if(length(idx)==0) next
  global_draws <- apply(agg_global_array[h_idx, idx, ], 2, mean, na.rm=TRUE)
  dom_draws <- apply(agg_domestic_array[h_idx, idx, ], 2, mean, na.rm=TRUE)
  period_summary$global_mean[i] <- median(global_draws, na.rm=TRUE)
  period_summary$domestic_mean[i] <- median(dom_draws, na.rm=TRUE)
  period_summary$P_global_gt_domestic[i] <- mean(global_draws > dom_draws, na.rm=TRUE)
}

print(period_summary)
write_xlsx(period_summary, "FEVD Data.xlsx")

# ================== Plots ==================
# FEVD shares over time
fevd_long <- fevd_df %>%
  pivot_longer(-date, names_to="type", values_to="value")
ggplot(fevd_long, aes(x=date, y=value, color=type)) +
  geom_line(size=0.8) +
  theme_minimal() +
  labs(title=paste("FEVD shares for", target_var, "(", choose_h, "-month horizon)"))

# Stacked decomposition
hist_decomp_df <- data.frame(date=rep(dates,3),
                             share=c(global_median_ts, dom_median_ts, 1-global_median_ts-dom_median_ts),
                             component=rep(c("Global","Domestic","Other"), each=length(dates)))
ggplot(hist_decomp_df, aes(x=date, y=share, fill=component)) +
  geom_area() +
  theme_minimal() +
  labs(title=paste("Median FEVD decomposition for", target_var, "(h=", choose_h, ")"))
