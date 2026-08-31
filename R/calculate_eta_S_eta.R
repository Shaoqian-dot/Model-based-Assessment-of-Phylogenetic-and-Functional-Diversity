# -----------------------------------------------------------------------------
# Calculate eta^T S eta for each observation
# -----------------------------------------------------------------------------
#
# Arguments:
#   abundance : matrix of observations x species
#   S         : species-by-species similarity matrix
#
# Returns:
#   A vector containing eta^T S eta for each observation.
#
# The species order in abundance must match the row/column order of S.

calculate_eta_S_eta <- function(abundance, S) {
  
  com <- abundance[, 1]
  # Check dimensions.
  if (
    ncol(abundance[, -1]) != nrow(S) ||
    ncol(abundance[, -1]) != ncol(S)
  ) {
    stop(
      "The number of species in abundance must match the dimensions of S."
    )
  }
  
  # Check species names if available.
  if (
    !is.null(colnames(abundance[, -1])) &&
    !is.null(rownames(S))
  ) {
    if (!all(colnames(abundance[, -1]) == rownames(S))) {
      stop(
        "Species order in abundance and S does not match."
      )
    }
  }
  
  # Calculate eta^T S eta for every row.
  abundanceData <- as.matrix(abundance[, -1])
  storage.mode(abundanceData) <- "numeric"
  result <- rowSums(
    (abundanceData %*% S) * abundanceData
  )
  DivModel <- data.frame(
    group = com,
    value = result
  )
  return(DivModel)
}