# Function to compute Rao's Q diversity and perform group comparisons
# abundance: matrix of species abundance (rows = samples, columns = species)
# DM_phy_func: distance matrix (functional or phylogenetic)
# randomization: indicator (0 = no randomization, 1 = apply randomization method)
# q: number of communities (groups)
safe_t_pvalue <- function(x, y, tol = 1e-12) {
  
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  tryCatch(
    t.test(
      x,
      y,
      paired = FALSE
    )$p.value,
    
    error = function(e) {
      
      # Stop for every other type of error
      if (!grepl(
        "data are essentially constant",
        conditionMessage(e),
        fixed = TRUE
      )) {
        stop(e)
      }
      
      x_constant <- max(abs(x - x[1])) <=
        tol * max(1, max(abs(x)))
      
      y_constant <- max(abs(y - y[1])) <=
        tol * max(1, max(abs(y)))
      
      # Stop if the error is not explained by both groups being constant
      if (!x_constant || !y_constant) {
        stop(e)
      }
      
      same_constant <- abs(x[1] - y[1]) <=
        tol * max(1, abs(x[1]), abs(y[1]))
      
      if (same_constant) {
        1
      } else {
        0
      }
    }
  )
}

getRaosQ <- function(abundance, DM_phy_func, use_randomization, q) {
  
  # Number of replicates per community
  r <- nrow(abundance) / q
  
  # Number of species
  m <- ncol(DM_phy_func)
  
  # Create community labels (1, 2, ..., q)
  com <- factor(rep(1:q, each = r))
  
  # Define groups to compare against the reference group (group 1)
  group_ids <- as.character(2:q)  # Comparisons: group 2–q vs group 1
  ref_group <- '1'
  
  # Case 1: No randomization (standard Rao's Q)
  if (use_randomization == 0){
    
    # Compute Rao's Q diversity for each sample
    # Output: vector of length r*q
    PD <- RaoQ(my.sample = abundance, DM = DM_phy_func)
    
    # Assign community labels as names
    names(PD) <- com
    
    # Initialize vector to store p-values
    p_value <- c()
    
    # Perform pairwise Wilcoxon rank-sum tests (unpaired)
    for (i in seq_along(group_ids)){
      p_value[i] <- wilcox.test(
        PD[names(PD) == ref_group],       # reference group (com 1)
        PD[names(PD) == group_ids[i]],    # comparison group
        paired = FALSE
      )$p.value 
    }
    
    # Case 2: With randomization
  } else if (use_randomization == 1) {
    
    # Compute Rao's Q (or equivalent metric) after randomization
    PD <- randomization(abundance = abundance, DM_phy_func = DM_phy_func)
    
    # Assign community labels
    names(PD) <- com
    
    # Initialize vector to store p-values
    p_value <- c()
    
    # Perform pairwise t-tests (unpaired)
    for (i in seq_along(group_ids)){
      x <- unname(PD[names(PD) == ref_group])
      y <- unname(PD[names(PD) == group_ids[i]])
      
      p_value[i] <- safe_t_pvalue(x, y)
      # p_value[i] <- t.test(
      #   PD[names(PD) == ref_group],       # reference group (com 1)
      #   PD[names(PD) == group_ids[i]],    # comparison group
      #   paired = FALSE
      # )$p.value
    }
  }
  
  # Assign names to the p-values for clarity
  names(p_value) <- paste0(2:q, "VS1")
  
  # Return vector of p-values for group comparisons
  return(p_value)
}