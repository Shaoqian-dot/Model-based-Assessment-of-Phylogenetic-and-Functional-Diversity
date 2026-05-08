# Create a function for finding NOPS pairs of species that cause no, small, or big diversity change if these pairs are swop
# We don't put phylogenetically same species in the the Index_Small_change or Index_Big_change, since if-else if-else if are mutually exclusive conditions
GroupingPairs <- function(DM_phy_func, m, quantile){
  Index_No_change <- c()
  Index_Small_change <- c()
  Index_Big_change <- c()
  #while(min(length(Index_No_change), length(Index_Small_change),
  #         length(Index_Big_change)) < (NOPS*2))
  {
    #shuffled <- sample(1:p)
    G1 <- seq(from = 1, to = m, by = 2)
    G2 <- seq(from = 2, to = m, by = 2)
    #Find pairs of species that cause no, small, or big diversity change if one or more pairs are swapped.
    dis_quant <- quantile(DM_phy_func[upper.tri(DM_phy_func)], quantile)
    for (Index_G1 in G1){
      for (Index_G2 in G2){
        {
          if (all.equal(unname(DM_phy_func[Index_G1, -Index_G1]), 
                        unname(DM_phy_func[Index_G2, -Index_G2]))==TRUE){
            while (!any(c(Index_G1, Index_G2) %in% Index_No_change)){
              Index_No_change <- c(Index_No_change, c(Index_G1, Index_G2))
            }
          }
          else if ((DM_phy_func[Index_G1, Index_G2] <= dis_quant) & (DM_phy_func[Index_G1, Index_G2] > 0)){
            while (!any(c(Index_G1, Index_G2) %in% Index_Small_change)){
              Index_Small_change <- c(Index_Small_change, c(Index_G1, Index_G2))
            }
          }
          else if (DM_phy_func[Index_G1, Index_G2] > dis_quant){
            while (!any(c(Index_G1, Index_G2) %in% Index_Big_change)){
              Index_Big_change <- c(Index_Big_change, c(Index_G1, Index_G2))
            }
          }
        }
      }
    }  
  }
  Index_change <- list(Index_No_change, Index_Small_change, Index_Big_change)
  assign("Index_change", Index_change, envir = .GlobalEnv)
  return(Index_change)
  # Find pairs of species that cause no, small, or big diversity change if one or more pairs are swapped.
}