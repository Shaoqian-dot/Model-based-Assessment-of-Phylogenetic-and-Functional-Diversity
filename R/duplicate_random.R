# Function to randomly duplicate elements in a vector
# x: input vector indicating indexes sampled from the full species pool in the real application data set
# m: desired final length after duplication
duplicate_random <- function(x, m) {
  # Number of elements to select before duplication
  m_less <- length(x)
  # Number of elements to randomly duplicate
  n <- m - m_less
  if (n > length(x)) stop("n cannot be greater than the length of the vector.")
  
  # Randomly select 'n' positions in 'x' to duplicate, and sort them
  pos <- sort(sample(seq_along(x), n, replace = FALSE))
  
  # Create a new vector to store the result with duplicates
  x_new <- x
  offset <- 0  # Offset accounts for the increasing length of x_new as we insert duplicates
  
  # Insert duplicates after the selected positions
  for (i in pos) {
    idx <- i + offset  # Adjust index to account for previous insertions
    x_new <- append(x_new, x_new[idx], after = idx) 
    offset <- offset + 1  # Increment offset after each insertion
  }
  
  # Return the new vector with random duplicates inserted
  return(x_new)
}