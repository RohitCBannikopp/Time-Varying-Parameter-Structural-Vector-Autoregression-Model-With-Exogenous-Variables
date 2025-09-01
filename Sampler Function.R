library(KFAS)
library(stochvol)
library(MASS)

run_tvp_svarx <- function(endo_clean, exo_clean, Z_array,
                          n_iter = 10000, burnin = 3000,
                          q_var = 1e-2,
                          progress_callback = NULL) {
  
  # Dimensions
  T <- nrow(endo_clean)
  n <- ncol(endo_clean)
  k_exo <- ncol(as.matrix(exo_clean))
  
  # Pre-allocate storage for draws
  theta_store <- vector("list", n_iter)
  h_store     <- vector("list", n_iter)
  A_store     <- vector("list", n_iter)
  eps_store   <- vector("list", n_iter)
  u_store     <- vector("list", n_iter)
  
  # Helper function to sample A_t states
  sample_A_state <- function(dep_var, regressors, sv_path, q_var_local) {
    sv_path[!is.finite(sv_path)] <- mean(sv_path, na.rm = TRUE)
    sv_path <- pmin(pmax(sv_path, -10), 10)
    H_array <- array(exp(sv_path), dim = c(1, 1, length(sv_path)))
    Z_array_A <- array(t(regressors), dim = c(1, ncol(regressors), length(dep_var)))
    
    model_A <- SSModel(dep_var ~ -1 + SSMcustom(
      Z = Z_array_A,
      T = diag(ncol(regressors)),
      Q = diag(q_var_local, ncol(regressors)),
      a1 = rep(0, ncol(regressors)),
      P1 = diag(1e4, ncol(regressors))
    ), H = H_array)
    
    state_draw <- simulateSSM(model_A, type = "states", conditional = TRUE)
    return(as.numeric(state_draw))
  }
  
  # ---------------- Gibbs Sampling ----------------
  for (iter in 1:n_iter) {
    
    theta_draws <- vector("list", n)
    
    # Sample coefficient paths (θ_t and γ_t)
    for (i in 1:n) {
      Z_i <- Z_array[, i, ]
      if (is.vector(Z_i)) Z_i <- matrix(Z_i, nrow = T)
      m_state <- ncol(Z_i)
      
      # Convert to 3D array for SSMcustom
      Zfull_arr <- array(NA, dim = c(1, m_state, T))
      for (tt in 1:T) Zfull_arr[1, , tt] <- Z_i[tt, ]
      
      # Q matrix: RW innovation variance
      # Note: In this implementation, exogenous coefficients share the same q_var
      # as endogenous coefficients.
      Q_mat <- diag(q_var, m_state)
      
      # Build state-space model
      model_i <- SSModel(
        endo_clean[, i] ~ -1 + SSMcustom(
          Z = Zfull_arr,
          T = diag(m_state),
          Q = Q_mat,
          a1 = rep(0, m_state),
          P1 = diag(1e6, m_state)  # diffuse prior
        ),
        H = array(1, dim = c(1, 1, T))
      )
      
      # Draw states using simulateSSM
      state_draw <- simulateSSM(model_i, type = "states", nsim = 1, conditional = TRUE)
      
      # Reshape to T × m_state
      state_mat <- matrix(state_draw[, 1, ], nrow = T, ncol = m_state, byrow = TRUE)
      theta_draws[[i]] <- list(alphahat = state_mat)
    }
    
    # Compute reduced-form residuals u_t
    u_mat <- matrix(NA, nrow = T, ncol = n)
    for (i in 1:n) {
      state_it <- theta_draws[[i]]$alphahat
      Z_i <- Z_array[, i, ]
      if (is.vector(Z_i)) Z_i <- matrix(Z_i, nrow = T)
      Y_hat <- rowSums(Z_i * state_it)
      u_mat[, i] <- endo_clean[, i] - Y_hat
    }
    
    # Sample stochastic volatilities
    h_t_mat <- matrix(NA, T, n)
    for (i in 1:n) {
      sv_post <- suppressWarnings(
        suppressMessages(
          svsample(u_mat[, i], draws = 1, burnin = 0,
                   priormu = c(-10, 1), priorphi = c(25, 5),
                   priorsigma = 0.1, quiet = TRUE)
        )
      )
      latent_data <- sv_post$latent
      if (is.list(latent_data)) latent_data <- latent_data[[1]]
      if (is.matrix(latent_data)) latent_data <- as.vector(latent_data)
      h_t_mat[, i] <- latent_data
    }
    
    # Sample A_t lower-triangular entries
    a_paths <- list()
    ij_index <- list()
    for (i in 2:n) {
      for (j in 1:(i - 1)) {
        path <- sample_A_state(u_mat[, i], cbind(u_mat[, j]), h_t_mat[, i], q_var)
        a_paths[[length(a_paths) + 1]] <- path
        ij_index[[length(ij_index) + 1]] <- c(i, j)
      }
    }
    
    # Build A_t list
    A_t_list <- vector("list", T)
    for (tt in 1:T) {
      A_now <- diag(1, n)
      for (k in seq_along(a_paths)) {
        idx <- ij_index[[k]]
        A_now[idx[1], idx[2]] <- -a_paths[[k]][tt]
      }
      A_t_list[[tt]] <- A_now
    }
    
    # Compute structural shocks ε_t
    eps_t <- sapply(1:T, function(tt) A_t_list[[tt]] %*% u_mat[tt, ])
    
    # Store draws
    theta_store[[iter]] <- theta_draws
    h_store[[iter]]     <- h_t_mat
    A_store[[iter]]     <- A_t_list
    eps_store[[iter]]   <- eps_t
    u_store[[iter]]     <- u_mat 
    
    # Progress callback
    if (iter %% 50 == 0 && !is.null(progress_callback)) {
      progress_callback(iter)
    }
  }
  
  return(list(
    theta = theta_store,
    h     = h_store,
    A     = A_store,
    eps   = eps_store,
    u     = u_store
  ))
}
