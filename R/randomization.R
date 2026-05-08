# Create a function to perform randomization and compute standardized effect size (SES)
randomization <- function(abundance, DM_phy_func){
  # Number of species
  m <- nrow(DM_phy_func)
  
  # Assign species names to abundance matrix columns
  colnames(abundance) <- colnames(DM_phy_func)
  
  # Generate null distribution of Rao's Q via tip randomization (999 permutations)
  null.output <- replicate(999, rand.RaoQ.fun(DM = DM_phy_func, my.sample = abundance), simplify = "matrix")
  
  # Compute observed Rao's Q using the original phylogeny
  # Ensure abundance columns match the phylogeny tip order
  obs_RaoQ <- RaoQ(
    my.sample  = abundance, 
    DM = DM_phy_func
  )
  
  # Compute standardized effect size (SES):
  # (observed - mean(null)) / sd(null)
  ses.all <- (obs_RaoQ -
                apply(null.output, MARGIN = 1, mean, na.rm = TRUE)) /
    apply(null.output, MARGIN = 1, sd, na.rm = TRUE)
  # Handle NaN values:
  # NaN occurs when sd(null) = 0 (no variation in null distribution)
  # In such cases, SES is set to 0 based on the assumption that
  # observed ≈ mean(null), hence no deviation from null expectation
  ses.all[is.nan(ses.all)] <- 0
  
  return(ses.all)
}

