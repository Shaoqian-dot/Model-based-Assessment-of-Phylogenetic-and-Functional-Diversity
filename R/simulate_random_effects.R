
# -----------------------------------------------------------------------------
# Simulate_random_effects
# -----------------------------------------------------------------------------
simulate_random_effects <- function(model) {
  # Conditional means of random effects
  ranef_model <- ranef(model) # condVar = TRUE is defaluted. 
  re_mean <- ranef_model$cond$id
  
  # Number of IDs and random-effect dimensions
  n <- nlevels(model$frame$id)
  p <- nlevels(model$frame$sp)
  
  if (ncol(re_mean) == p){
    # Conditional variance-covariance matrices
    re_var <- attr(re_mean, "condVar")
  } else if(ncol(re_mean) == (p + 1)) {
    # Conditional variance-covariance matrices
    re_var <- attr(re_mean, "condVar")
    re_var_Intercept <- re_var$`1`
    re_var_sp <- re_var$`2`
  }
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
    
    # Generate one random-effect vector
    if (length(mu_j) == p){
      
      # Conditional covariance for ID j
      Sigma_j <- re_var[, , j]
      Sigma_j[is.na(Sigma_j)] <- 0
      
      RR[j, ] <- MASS::mvrnorm(
        n = 1,
        #mu = mu_j,
        mu = rep(0, times = length(mu_j)),
        Sigma = Sigma_j
      )
    } else if (length(mu_j) == (p + 1)){
      # Conditional covariance for ID j
      Sigma_j <- re_var_sp[, , j]
      Sigma_j[is.na(Sigma_j)] <- 0
      
      RR[j, ] <- MASS::mvrnorm(
        n = 1,
        #mu = mu_j,
        mu = rep(0, times = p),
        Sigma = Sigma_j
      )  + rnorm(1, mean = 0, sd = sqrt(re_var_Intercept[, , j]))
      
    }
    # Keep random-effect names
    colnames(RR) <- paste0('sp', levels(model$frame$sp))
  }
  return(RR)
}
