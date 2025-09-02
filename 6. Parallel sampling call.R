# Load required libraries
library(future.apply)   # for parallel execution
library(progressr)      # for progress tracking
library(parallel)       # for detecting CPU cores
library(KFAS)           # state space modeling
library(stochvol)       # stochastic volatility models
library(MASS)           # matrix utilities

# --- Setup parallel execution ---
# Use all but one available core
plan(multisession, workers = detectCores() - 1)

# Choose a simple progress handler (prints a progress bar in console)
handlers("txtprogressbar")

# Set seed for reproducibility (important in Bayesian MCMC)
set.seed(12345)

# --- Run chains in parallel with progress monitoring ---
with_progress({
  
  # Total updates = number of chains * (iterations / reporting frequency)
  # Example: 3 chains * (10000 iterations / 50 update interval) = 600 updates
  total_updates <- 3 * (10000 / 50)  
  p <- progressor(along = seq_len(total_updates))  # initialize progress tracker
  
  # Run 3 parallel chains using future_lapply
  chains1 <- future_lapply(1:3, function(chain_id) {
    
    # Call user-defined function run_tvp_svarx()
    run_tvp_svarx(
      endo_clean, exo_clean, Z_array, models,   # input data and model specs
      n_iter = 10000,                           # number of MCMC iterations
      burnin = 3000,                            # burn-in samples to discard
      q_var = 1e-2,                             # state noise variance parameter
      progress_callback = function(iter) {      # callback for progress updates
        p(sprintf("Chain %d: Iter %d", chain_id, iter))
      }
    )
    
  }, future.seed = TRUE)   # ensure reproducibility across parallel workers
})

# --- Reset to sequential execution ---
plan(sequential)
