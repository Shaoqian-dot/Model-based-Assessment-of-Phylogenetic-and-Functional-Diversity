# Assumptions:
# 1. VC_phy_func is a square matrix
# 2. VC_phy_func is symmetric (required for real eigenvalues and orthogonal eigenvectors)
# 3. No missing values (NA/NaN) are present in the matrix

# Function to perform spectral (eigen) decomposition of a matrix
spectral_decomp <- function(VC_phy_func){
  
  # Step 1: Compute the eigen decomposition of the input matrix
  # VC_phy_func is assumed to be a square matrix (e.g., covariance or similarity matrix)
  # The decomposition returns:
  #   - eigenvalues (spec_decomp$values)
  #   - eigenvectors (spec_decomp$vectors)
  spec_decomp <- eigen(VC_phy_func)
  
  # Step 2: Extract the matrix of eigenvectors
  # Columns of P are eigenvectors of VC_phy_func
  P <- spec_decomp$vectors
  
  # Step 3: Construct a diagonal matrix of eigenvalues
  # Eigenvalues are placed along the diagonal of matrix D
  D <- diag(spec_decomp$values) 
  
  # Step 4: Return the decomposition components
  # The original matrix can be reconstructed as:
  #   VC_phy_func = P %*% D %*% t(P)   (if symmetric)
  return(list(P = P, D = D))
}