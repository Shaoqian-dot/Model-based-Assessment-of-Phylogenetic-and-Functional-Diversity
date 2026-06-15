# Function to simulate species abundance data based on a functional distance matrix
# DM_phy_func: functional distance matrix of species
# NOPS: number of pair swaps (used in mean_matrix)
# q: number of communities
# p: number of species
# alpha, beta: parameters controlling mean abundance
# r: number of replicates per community
# quantile: threshold for defining small/large diversity changes
# Corr, P, Eigen: parameters for mean_matrix function
getData <- function (DM_phy_func, NOPS, q, 
                     alpha, beta, r, quantile,
                     Corr, P, Eigen, Distribution, p){
  
  # Determine which species pairs will have diversity changes (small/no/big)
  # Returns a matrix/list of indices for changing species pairs
  Index_change <- GroupingPairs(DM_phy_func = DM_phy_func, quantile = quantile, p = p)

  # Generate the mean abundance matrix for each community and species
  # Uses parameters alpha, beta, q, m, r, NOPS, Corr, P, Eigen
  Mean_matrix_full <- mean_matrix(Index_change = Index_change, beta = beta,
                                  alpha = alpha, q = q, r = r, NOPS = NOPS,
                                  Corr = Corr, P = P, Eigen = Eigen, Distribution = Distribution, p = p)
  # Simulate species abundance using a Bernoulli/binomial process
  # Each element is 0 or 1, probability given by corresponding element in Mean_matrix_full
  if (Distribution == 'binomial') {
    
    abundance <- matrix(
      rbinom(length(Mean_matrix_full), 1, as.vector(Mean_matrix_full)),
      ncol = p
    )
    
  } else if (Distribution == 'poisson') {
    
    abundance <- matrix(
      rpois(length(Mean_matrix_full), as.vector(Mean_matrix_full)),
      ncol = p
    )
    
  }
  
  
  # Name the species columns as "species.1", "species.2", ...
  # Note that "species.1", "species.2", ... correspond to colnames(DM_phy_func)/rownames(DM_phy_func), respectively. 
  # We use species 1/sp 1,...,species p/sp p here for convenience.
  colnames(abundance) <- paste("species", 1:p, sep = ".")
  
  # Create community factor (1:q repeated r times each)
  com <- factor(rep(1:q, each=r))
  
  # Combine abundance matrix and community labels into a data frame
  yX <- data.frame(abundance, com)
  
  # Reshape data frame from wide to long format:
  #   y: abundance value
  #   sp: species
  yX <- reshape(yX, direction = "long", varying = colnames(abundance),
                v.names = "y", timevar = "sp")
  
  # Convert species column to factor
  yX$sp <- factor(yX$sp)
  
  # Return a list with reshaped long-format data frame and the abundance matrix
  return(list(yX = yX, abundance = abundance))
}