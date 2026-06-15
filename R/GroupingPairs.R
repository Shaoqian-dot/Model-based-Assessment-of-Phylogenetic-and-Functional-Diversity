# The idea here is to separate 
GroupingPairs <- function(DM_phy_func, quantile, p){
  # Number of species
  m <- p
  
  # Initialize lists (use list for efficiency, convert later if needed)
  Index_No_change   <- integer(0)
  Index_Small_change <- integer(0)
  Index_Big_change   <- integer(0)
  
  # Define groups
  G1 <- seq(1, m, by = 2)
  G2 <- seq(2, m, by = 2)
  
  # Compute distance threshold
  dis_quant <- quantile(DM_phy_func[upper.tri(DM_phy_func)], quantile)
  
  for (i in G1){
    for (j in G2){
      
      # Case 1: identical phylogenetic profile
      if (isTRUE(all.equal(
        unname(DM_phy_func[i, -i]),
        unname(DM_phy_func[j, -j])
      ))) {
        
        if (!any(c(i, j) %in% Index_No_change)) {
          Index_No_change <- c(Index_No_change, i, j)
        }
        
        # Case 2: small change
      } else if (DM_phy_func[i, j] <= dis_quant) { # Pairs in small change group exclude those pairs in no change group since if-else if are excluded. 
        
        if (!any(c(i, j) %in% Index_Small_change)) {
          Index_Small_change <- c(Index_Small_change, i, j)
        }
        
        # Case 3: big change
      } else if (DM_phy_func[i, j] > dis_quant) { # Pairs in large change group exclude those pairs in no change group since if-else if are excluded. 
        
        if (!any(c(i, j) %in% Index_Big_change)) {
          Index_Big_change <- c(Index_Big_change, i, j)
        }
      }
    }
  }
  
  # Return result (no global assignment!)
  return(list(
    No_change   = Index_No_change,
    Small_change = Index_Small_change,
    Big_change   = Index_Big_change
  ))
}
