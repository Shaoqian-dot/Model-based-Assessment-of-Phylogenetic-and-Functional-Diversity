# Function to compute Rao's Q after randomly shuffling species labels
# on the phylogenetic (or trait-based) distance matrix
rand.RaoQ.fun <- function(DM, my.sample){ 
  
  # Extract species names from the distance matrix (column names)
  species <- colnames(DM)
  
  # Generate a random permutation of species labels (without replacement)
  species_permutation <- sample(species, size = length(species), replace = FALSE)
  
  # Shuffle both rows and columns of labels of the distance matrix according to the permutation
  # This preserves the distance structure but randomizes species identities.
  # Note: Permuting species labels on the distance matrix is mathematically equivalent 
  # to permuting labels on the underlying trait or phylogenetic data and recomputing 
  # the distance matrix, provided that the permutation is applied consistently to 
  # both rows and columns.
  DM_shuffle <- DM
  colnames(DM_shuffle) <- species_permutation
  rownames(DM_shuffle) <- species_permutation
  
  # Reorder the community matrix to match the shuffled distance matrix
  # and compute Rao's quadratic entropy
  RaoQ(
    my.sample = my.sample[, colnames(DM_shuffle)], 
    DM = DM_shuffle
  )
}