U <- function(r, q, P, Eigen){
  # Number of species
  m <- ncol(P)
  
  # Initialize U matrix of size (r*q) x m
  U <- matrix(NA, r * q, m)
  
  
  # Select which pair of eigenvectors to use
  if (Eigen == 1){
    Lambda <- P[, c(1, 2)]   # Use the first two eigenvectors
    #rep <- Model_fix_J_rr$obj$env$report()
    #Lambda <- rep$fact_load[[1]]
    #Lambda <- matrix(c(0.34641103,  0.40831772,  0.43208410,  0.58132042, -0.05318892,  
    #                   0.00000000,  0.53786800,  0.27269484,  0.48162548, -0.19023844), m, 2)
  } else if(Eigen == 2) {
    
    Lambda <- P[, c(m - 1, m)]  # Use the last two eigenvectors
  
    } else if (Eigen == 3) {
      
    mu_hat <- c(-0.21, 0.42)
    Sigma_hat <-matrix(c(0.61, 0.06, 
                          0.06, 0.99), 2, 2)
    Lambda <- mvrnorm(m, mu_hat, Sigma_hat)
  } else if (Eigen == 4) {
    
    mu_hat <- c(-0.21, 0.42)
    Sigma_hat <- 1.5 * matrix(c(0.61, 0.06, 
                          0.06, 0.99), 2, 2)
    Lambda <- mvrnorm(m, mu_hat, Sigma_hat)
  }
  
  # Loop over all rows (r*q) and species (m)
  for (i in 1 : (r * q)){
    z <- mvrnorm(1, rep(0, 2), diag(c(1, 1)))
    for (j in 1 : m){
      # Generate a random 2D vector from standard normal distribution
      # Multiply it by the j-th row of Lambda to get a random noise
      # Store the scalar in U[i, j]
      U[i, j] <- z %*% Lambda[j, ]
    }
  }
  
  # Return the simulated contrast matrix
  return(U)
}