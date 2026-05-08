# Function to construct a functional distance matrix for a subset of species
# with random duplication to avoid identical mean vectors in later pair-swap steps
# trait: data frame or matrix of species traits
# m: desired number of species after duplication
# m_less: number of unique species to sample before duplication
getFuncMatrix <- function(trait, m, m_less){
  
  # Compute the Gower distance matrix based on species traits
  DM_trait <- as.matrix(gowdis(trait))
  
  # Total number of species available
  n_species_total <- nrow(trait)
  
  # Randomly select m_less unique species (without replacement)
  index <- sample(1 : n_species_total, m_less, replace = FALSE)
  
  # Randomly duplicate (m - m_less) species among the selected ones.
  # Duplicated species are inserted immediately after the originals.
  # Purpose:
  #   - Avoid identical mean vectors across communities (e.g., com 3 vs com 4)
  #     caused by pair-swap operations.
  index_rep <- duplicate_random(x = index, m = m)
  
  # Extract the functional distance submatrix for the selected (and duplicated) species
  DM_phy_func <- DM_trait[index_rep, index_rep]
  
  # Return the resulting functional distance matrix
  return(DM_phy_func)
}