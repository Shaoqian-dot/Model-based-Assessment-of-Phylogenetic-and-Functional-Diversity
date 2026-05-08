nameContrast <- function (Data, val_num.eig){
  m <- length(levels(Data$sp))
  # Number of contrasts for a factor with m levels is (m - 1)
  n_contr <- m - 1
  
  {
    # If the number of eigenvectors is less than the number of contrasts
    if (val_num.eig < n_contr) {
      
      # Assign column names to the contrast matrix:
      # First part: "C1", "C2", ..., up to the number of eigenvectors
      # Remaining part: "C_S..." indicating supplementary contrasts
      colnames(contrasts(Data$sp)) <- c(
        paste0("C", seq_len(val_num.eig)),
        paste0("C_S", seq(val_num.eig + 1, n_contr))
      )
      
    } else {
      
      # If enough eigenvalues are available, name all contrasts as "C1", ..., "C(n_contr)"
      colnames(contrasts(Data$sp)) <- paste0("C", seq_len(n_contr))
    } 
  }
  
  # Return the modified Data object with updated contrast names
  return(Data)
}
