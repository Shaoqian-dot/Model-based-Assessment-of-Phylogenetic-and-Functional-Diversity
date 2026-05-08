# vec: the original vector (length must be even)
# n_sample: number of pairs to sample

sample_pairs <- function(vec, n_sample){
  # Indices of the first element in each pair (1, 3, 5, ...)
  odd_numbers <- seq(1, length(vec), by = 2)                
  
  # Total number of pairs
  n_pairs <- length(vec) / 2 
  
  # Randomly select indices corresponding to pairs
  selected_pair_idx <- sample(odd_numbers, n_sample)
  
  # Get the indices of both elements in each selected pair
  selected_elements_idx <- as.vector(
    sapply(selected_pair_idx, function(i) c(i, i + 1))
  )
  
  # Return the selected elements
  return(vec[selected_elements_idx])
}