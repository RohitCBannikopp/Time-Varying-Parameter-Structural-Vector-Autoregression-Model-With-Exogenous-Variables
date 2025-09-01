library(dplyr)
library(ggplot2)
library(coda)
library(viridis)
library(cowplot)

# Function to compute time-varying R-hat
compute_rhat_df <- function(chains, burnin = 3000, eq_names = NULL) {
  
  n_chains <- length(chains)
  n_iter   <- length(chains[[1]]$theta)
  n_eq     <- length(chains[[1]]$theta[[1]])
  
  if (is.null(eq_names)) {
    eq_names <- paste("Equation", seq_len(n_eq))
  }
  
  rhat_df <- bind_rows(lapply(seq_len(n_eq), function(eq) {
    alphahat_example <- chains[[1]]$theta[[1]][[eq]]$alphahat
    T       <- nrow(alphahat_example)
    n_coef  <- ncol(alphahat_example)
    
    bind_rows(lapply(seq_len(n_coef), function(ci) {
      rhat_ts <- sapply(seq_len(T), function(t) {
        draws <- lapply(seq_len(n_chains), function(ch) {
          sapply((burnin + 1):n_iter, function(i) {
            chains[[ch]]$theta[[i]][[eq]]$alphahat[t, ci]
          })
        })
        mcmc_list <- as.mcmc.list(lapply(draws, as.mcmc))
        gelman.diag(mcmc_list, autoburnin = FALSE)$psrf[1]
      })
      
      tibble(
        equation    = eq_names[eq],
        coef_index  = ci,
        coefficient = if (ci == 1) "Intercept" else paste0("coef.", ci - 1),
        time        = seq_len(T),
        Rhat        = rhat_ts
      )
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
