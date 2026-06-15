# Create a function for producing a mean matrix
# U() must be defined externally and return a matrix of dimension (q*r) x m
mean_matrix <- function (Index_change, beta, alpha, 
                         q, r, NOPS, Corr, P, Eigen, Distribution){
  # Number of species
  m <- nrow(P)
  
  # Construct (q-1) rows, each with m elements, alternating between -beta and beta
  Means_logit_beta <- rep(rep(c(-beta, beta), length.out = m), times = q - 1)
  
  # Construct the last row (community 5), each with m elements,
  # alternating between (alpha - beta) and (alpha + beta)
  Means_logit_alpha <- rep(c(alpha - beta, alpha + beta), length.out = m)
  
  # Concatenate all rows into a single vector
  Means_logit <- c(Means_logit_beta, Means_logit_alpha)
  
  # Reshape into a q x m matrix (row-wise filling)
  Mean_matrix_logit <- matrix(Means_logit, c(q, m), byrow = TRUE)
  
  # Loop over groups (excluding the first and last group)
  for (i in 2 : (q - 1)){
    
    # Number of available index pairs for swapping in group (i-1)
    available_pairs <- length(Index_change[[i-1]]) / 2
    
    # If no valid pairs exist, print a message
    if(available_pairs == 0){
      message(sprintf("No satisfied pair in group %d", i-1))
    } else {
      # Number of pairs to sample (cannot exceed available pairs)
      n_sample <- min(NOPS, available_pairs)
      
      # Sample index pairs from the provided index set
      Pairs <- sample_pairs(vec = Index_change[[i-1]], n_sample = n_sample)
      
      # Flip the sign of selected positions in the mean matrix
      Mean_matrix_logit[i, Pairs] <- -Mean_matrix_logit[i, Pairs]
    }
  }
  
  # If correlation structure is required
  if (Corr == 1){
    
    # Generate random effect matrix U_M of dimension (q*r) x m
    U_M <- U(p = m, r = r, q = q, P = P, Eigen = Eigen)
    
    # Expand mean matrix to (q*r) x m and add random effects
    # Then apply inverse logit transformation
    Mean_matrix_logit_cor <- Mean_matrix_logit[rep(1 : q, each = r), ] + U_M
    
  } else {
    # No correlation: only expand and apply inverse link function
    Mean_matrix_logit_cor <- Mean_matrix_logit[rep(1 : q, each = r), ]
  }
  
  if (Distribution == 'binomial') {
    Mean_matrix_full <- plogis(Mean_matrix_logit_cor)
  } else if (Distribution == 'poisson'){
    Mean_matrix_full <- exp(Mean_matrix_logit_cor)
  }
  assign('Mean_matrix_full', Mean_matrix_full, envir =  .GlobalEnv)
  # Return final mean matrix
  return(Mean_matrix_full)
}