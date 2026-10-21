U <- function(p, r, q, V_J = NA, Eigen, sp_Sim = NA, FactorLoading = NA){
  
  # Initialize U matrix of size (r*q) x p
  U <- matrix(NA, r * q, p)
  colnames(U) <- sp_Sim
  
  # Select which pair of eigenvectors to use
  if (Eigen == 1){
    
    Lambda <- V_J[, c(1, 2)]   # Use the first two eigenvectors
    
      } else if(Eigen == 2) {
    
    Lambda <- V_J[, (ncol(V_J)-1):ncol(V_J)]  # Use the last two eigenvectors
  
  } else if (Eigen == 3) {
    Lambda <- FactorLoading[sp_Sim, ]
    } 

  # Loop over all rows (r*q) and species (p)
  for (i in 1 : (r * q)){
    z <- mvrnorm(1, rep(0, 2), diag(c(1, 1)))
    for (j in 1 : p){
      # Generate a random 2D vector from standard normal distribution
      # Multiply it by the j-th row of Lambda to get a random noise
      # Store the scalar in U[i, j]
      U[i, j] <- z %*% Lambda[j, ]
    }
  }
  
  # Return the simulated contrast matrix
  return(U)
}