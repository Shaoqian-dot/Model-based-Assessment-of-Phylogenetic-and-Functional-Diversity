
# -----------------------------------------------------------------------------
# Simulate_random_effects
# -----------------------------------------------------------------------------
simulate_random_effects <- function(model) {
  
  # Conditional means of random effects
  ranef_model <- ranef(model, condVar = TRUE)
  re_mean <- ranef_model$cond$id
  
  # Conditional variance-covariance matrices
  re_var <- attr(re_mean, "condVar")
  
  # Number of IDs and random-effect dimensions
  n <- nrow(re_mean)
  p <- ncol(re_mean)
  
  # Matrix to store simulated random effects
  RR <- matrix(
    NA_real_,
    nrow = n,
    ncol = p
  )
  
  # Generate one multivariate-normal draw for each ID
  for (j in seq_len(n)) {
    
    # Conditional mean for ID j
    mu_j <- as.numeric(re_mean[j, ])
    
    # Conditional covariance for ID j
    Sigma_j <- re_var[, , j]
    Sigma_j[is.na(Sigma_j)] <- 0
    
    # Generate one random-effect vector
    RR[j, ] <- MASS::mvrnorm(
      n = 1,
      mu = mu_j,
      Sigma = Sigma_j
    )
  }
  
  # Keep random-effect names
  colnames(RR) <- colnames(re_mean)
  
  return(RR)
}
