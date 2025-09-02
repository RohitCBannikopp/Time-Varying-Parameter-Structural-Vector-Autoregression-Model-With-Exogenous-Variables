library(dplyr)
library(ggplot2)
library(coda)
library(viridis)
library(cowplot)

# Function to compute time-varying R-hat
compute_rhat_df <- function(chains, burnin = 3000, eq_names = NULL) {
  
  # Number of independent MCMC chains
  n_chains <- length(chains)
  # Number of iterations in each chain
  n_iter   <- length(chains[[1]]$theta)
  # Number of equations in the model
  n_eq     <- length(chains[[1]]$theta[[1]])
  
  # If no equation names are provided, create generic labels
  if (is.null(eq_names)) {
    eq_names <- paste("Equation", seq_len(n_eq))
  }
  
  # Loop over each equation
  rhat_df <- bind_rows(lapply(seq_len(n_eq), function(eq) {
    
    # Use the first chain, first iteration, current equation
    # just to extract the structure of alphahat (matrix of coefficients over time)
    alphahat_example <- chains[[1]]$theta[[1]][[eq]]$alphahat
    
    # T = number of time periods
    T       <- nrow(alphahat_example)
    # n_coef = number of coefficients for this equation
    n_coef  <- ncol(alphahat_example)
    
    # Loop over coefficients
    bind_rows(lapply(seq_len(n_coef), function(ci) {
      
      # For each coefficient, compute R-hat over time
      rhat_ts <- sapply(seq_len(T), function(t) {
        
        # Collect posterior draws across chains for this time and coefficient
        draws <- lapply(seq_len(n_chains), function(ch) {
          sapply((burnin + 1):n_iter, function(i) {
            chains[[ch]]$theta[[i]][[eq]]$alphahat[t, ci]
          })
        })
        
        # Convert draws to coda format (mcmc.list) for R-hat computation
        mcmc_list <- as.mcmc.list(lapply(draws, as.mcmc))
        
        # Compute Gelman-Rubin statistic (R-hat)
        gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[1]
      })
      
      # Return a tidy tibble with results
      tibble(
        equation    = eq_names[eq],                    # name of the equation
        coef_index  = ci,                              # coefficient index
        coefficient = if (ci == 1) "Intercept"         # label intercept
                       else paste0("coef.", ci - 1),   # label slope coefficients
        time        = seq_len(T),                      # time index
        Rhat        = rhat_ts                          # computed R-hat values
      ) %>%
        # Ensure coefficients have a consistent factor order
        mutate(coefficient = factor(
          coefficient,
          levels = c("Intercept", paste0("coef.", 1:(n_coef - 1)))
        ))
    }))
  }))
  
  return(rhat_df)
}

# ==== Example Usage ====
eq_labels <- c("N_50", "Real Rate", "Term Spread", "DP Ratio", "Relative bill rate")
rhat_df <- compute_rhat_df(chains1, burnin = 3000, eq_names = eq_labels)

# Define colors
all_coefs <- unique(rhat_df$coefficient)
coef_colors <- setNames(hcl.colors(length(all_coefs), "Dark3"), all_coefs)

# Plot R-hat per equation
plots <- lapply(unique(rhat_df$equation), function(eq_name) {
  df_eq <- filter(rhat_df, equation == eq_name, !is.na(Rhat))
  ggplot(df_eq, aes(x = time, y = Rhat, color = coefficient)) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_hline(yintercept = 0.9999, linetype = "dotted", color = "gray40") +
    geom_hline(yintercept = 1.0001, linetype = "dashed", color = "red") +
    scale_color_manual(values = coef_colors) +
    labs(title = paste0("R-hat Over Time for ", eq_name),
         x = "Time", y = "Gelman–Rubin R-hat") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
})

# Extract and plot legend
legend_plot <- ggplot(rhat_df, aes(x = time, y = Rhat, color = coefficient)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = coef_colors) +
  labs(color = "Coefficient") +
  theme_minimal(base_size = 13) +
  guides(color = guide_legend(override.aes = list(linewidth = 2)))
legend_only <- cowplot::get_legend(legend_plot)

cowplot::plot_grid(legend_only)

# Show plots
for (p in plots) print(p)

# Summarise R-hat values by equation
rhat_summary_eq <- rhat_df %>%
  filter(!is.na(Rhat)) %>%
  group_by(equation) %>%
  summarise(
    min_R  = round(min(Rhat, na.rm = TRUE), 6),
    max_R  = round(max(Rhat, na.rm = TRUE), 6),
    mean_R = round(mean(Rhat, na.rm = TRUE), 6),
    sd_R   = round(sd(Rhat, na.rm = TRUE), 6),
    .groups = "drop"
  )
print(rhat_summary_eq)
