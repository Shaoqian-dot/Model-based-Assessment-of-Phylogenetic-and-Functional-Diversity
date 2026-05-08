# Function to compute Rao's Quadratic Entropy (Rao's Q)
# Note that species names in my.sample and DM are consistent.
RaoQ <- function(my.sample, DM){  # DM: species-by-species distance matrix
  
  # Step 1: Convert abundance data to relative abundance (row-wise normalization)
  # Each row represents a site; values are scaled to sum to 1
  relative.my.sample <- my.sample / rowSums(my.sample)
  
  # Step 2: Handle rows with zero total abundance (division by zero → NaN)
  # Replace NaN values with 0, meaning no contribution to diversity
  relative.my.sample[is.na(relative.my.sample)] <- 0
  
  # Step 3: Store relative abundance matrix
  # P is an n × m matrix (n sites, m species)
  P <- relative.my.sample
  
  # Step 4: Compute Rao's Q using a fully vectorized formulation
  # Rao's Q for each site i is:
  #   Q_i = 0.5 * p_i^T * DM * p_i
  # This is efficiently computed for all sites as:
  #   0.5 * rowSums((P %*% DM) * P)
  RaoQ.diversity <- 0.5 * rowSums((P %*% DM) * P)
  
  # Step 5: Return a vector of Rao's Q values (one per site)
  return(RaoQ.diversity)
}
