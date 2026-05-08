PD_FDSimulation <- function (alpha, beta, q, m, m_less, NOPS, r, diversity,
                             Methods, Methods_para, tree, trait, quantile,
                             Corr, n_Model, Eigen, Distribution){ # num.eig can be replaced with `method`.
  Model_label <- c("fix", "fix_rr",
                   "p_J", "p_J_rr", 
                   "p_fix_J", "p_fix_J_rr", 
                   "fix_J", "fix_J_rr")
  Method_label <- c("EVs 1–2", "λ > mean")
  # Initialization of output list
  p_summary <- list()
  
  # 3D array to store p-values: (community × model × method)
  p_Model_Array <- array(NA, dim = c(q - 1, n_Model, length(Methods)))
  
  # Assign dimension names for clarity
  dimnames(p_Model_Array) <- list(
    row    = paste0("com", 2:q),         # communities (excluding baseline)
    col    = Model_label, # model indices
    face   = Method_label  # eigenvalue selection methods
  )
  
  # Construct distance matrix depending on diversity type
  if (diversity == 1) {
    # Phylogenetic distance matrix
    DM_phy_func <- getPhyloMatrix(tree = tree, m = m)
  } else if (diversity == 2) {
    # Functional trait distance matrix
    DM_phy_func <- getFuncMatrix(trait = trait, m = m, m_less = m_less)
  }
  
  # Convert distance matrix to similarity (variance-covariance-like) matrix
  VC_phy_func <- max(DM_phy_func) - DM_phy_func
  
  # Spectral decomposition (eigenvectors and eigenvalues)
  spec <- spectral_decomp(VC_phy_func = VC_phy_func)
  P <- spec$P  # eigenvectors
  D <- spec$D  # eigenvalues
  # Determine number of eigenvectors to use under different methods
  num.eig <- get_num.eig(Methods = Methods, Methods_para = Methods_para, D = D)
  
  # Generate simulated data (contains randomness)
  Data <- getData(DM_phy_func = DM_phy_func, NOPS = NOPS,
                  q = q, alpha = alpha, 
                  beta = beta, r = r, quantile = quantile,
                  Corr = Corr, P = P, Eigen = Eigen, Distribution = Distribution)

  yX <- Data$yX                 # response + predictors
  abundance <- Data$abundance   # species abundance matrix
  
  # ===== Precompute transformed similarity matrix (P_J) for Models 5–8 =====
  I_m <- diag(m) 
  J <- I_m - 1/m * matrix(1, m, m)   # centering matrix
  S_J <- t(J) %*% VC_phy_func %*% J  # centered similarity matrix
  P_J <- spectral_decomp(VC_phy_func = S_J)$P  # eigenvectors of transformed matrix
  
  # Matrix to track model failures (NA models due to singular contrast matix)
  NA_contrast <- matrix(NA, nrow = length(Methods), ncol = n_Model)
  dimnames(NA_contrast) <- list(
    row    = Method_label,
    column = Model_label
  )
  
  # ===== Loop over eigenvalue selection methods =====
  Distribution_negBino <- nbinom2()
  for (k in 1 : length(Methods)) {
    val_num.eig <- num.eig[k]  # number of eigenvectors used
    
    # ===== Fit 10 different model variants =====
    Model_fix <- fit_model(type = 'contrast', matrix_type = 'P', rr = FALSE,
                          yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    assign(' Model_fix',  Model_fix, env = .GlobalEnv)
    Model_fix_rr <- fit_model(type = 'contrast', matrix_type = 'P', rr = TRUE,
                             yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    
    # Model_p <- fit_model(type = 'p', matrix_type = 'P', rr = FALSE,
    #                     yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution)
    # 
    # Model_p_rr <- fit_model(type = 'p', matrix_type = 'P', rr = TRUE,
    #                        yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution)

    Model_p_J <- fit_model(type = 'p', matrix_type = 'P_J', rr = FALSE,
                          yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    
    Model_p_J_rr <- fit_model(type = 'p', matrix_type = 'P_J', rr = TRUE,
                             yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    
    Model_p_fix_J <- fit_model(type = 'p_fix', matrix_type = 'P_J', rr = FALSE,
                           yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    
    Model_p_fix_J_rr <- fit_model(type = 'p_fix', matrix_type = 'P_J', rr = TRUE,
                              yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    
    Model_fix_J <- fit_model(type = 'contrast', matrix_type = 'P_J', rr = FALSE,
                            yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    
    Model_fix_J_rr <- fit_model(type = 'contrast', matrix_type = 'P_J', rr = TRUE,
                               yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, Distribution = Distribution_negBino)
    # Function to detect model being NA
    is_failed_model <- function(x) {
      length(x) == 1 && is.na(x)
    }
    
    # Combine models into a list (1 × 10)
    models_list <- list(
      Model_fix, Model_fix_rr,
      Model_p_J, Model_p_J_rr, 
      Model_p_fix_J, Model_p_fix_J_rr, 
      Model_fix_J, Model_fix_J_rr
    )
    
    # NA_contrast: Record NA models caused by the singular contrast 
    NA_contrast[k, ] <- as.numeric(sapply(models_list, is_failed_model))
    
    # ===== Compute chi-square statistics =====
    chi_square <- matrix(NA, q - 1, n_Model)
    
    dimnames(chi_square) <- list(
      row    = paste0("com", 2:q),
      col    = Model_label
    )
    
    for(j in 2 : q){
      
      # Labels for extracting coefficients
      label_Model_1278 <- paste0('spC', (1 : val_num.eig), ':com', j)
      label_Model_3456 <- paste0('com', j, ':p', 1 : val_num.eig) 
      
      # Compute chi-square statistics for each model
      chi_square_1 <- getChisquare_fixed (Model = Model_fix,          label = label_Model_1278)
      chi_square_2 <- getChisquare_fixed (Model = Model_fix_rr,       label = label_Model_1278)
      #chi_square_3 <- getChisquare_fixed (Model = Model_p,            label = label_Model_345678)
      #chi_square_4 <- getChisquare_fixed (Model = Model_p_rr,         label = label_Model_345678)
      chi_square_3 <- getChisquare_fixed (Model = Model_p_J,          label = label_Model_3456)
      chi_square_4 <- getChisquare_fixed (Model = Model_p_J_rr,       label = label_Model_3456)
      chi_square_5 <- getChisquare_fixed (Model = Model_p_fix_J,      label = label_Model_3456)
      chi_square_6 <- getChisquare_fixed (Model = Model_p_fix_J_rr,   label = label_Model_3456)
      chi_square_7 <- getChisquare_fixed (Model = Model_fix_J,        label = label_Model_1278)
      chi_square_8 <- getChisquare_fixed (Model = Model_fix_J_rr,     label = label_Model_1278)
      
      # Store results in matrix
      chi_square[j - 1, ] <- c(
        chi_square_1, chi_square_2, 
        chi_square_3, chi_square_4, 
        chi_square_5, chi_square_6,
        chi_square_7, chi_square_8
      )
    }
    
    # Convert chi-square to p-values
    p_Model_Matrix <- 1 - pchisq(chi_square, df = val_num.eig)
    
    # Store results for this method
    p_Model_Array [, , k] <- p_Model_Matrix
    
    # Clean up large objects to save memory
    rm(Model_fix, Model_fix_rr,
       Model_p_J, Model_p_J_rr,
       Model_p_fix_J, Model_p_fix_J_rr,
       Model_fix_J, Model_fix_J_rr
    )
    gc()
  }
  
  # ===== Rao's Q diversity test =====
  p_RaosQ <- getRaosQ(abundance = abundance, DM_phy_func = DM_phy_func,
                      use_randomization = 0, q = q)
  
  # Randomized Rao's Q
  p_RaosQ_rand <- getRaosQ(abundance = abundance, DM_phy_func = DM_phy_func,
                           use_randomization = 1, q = q)
  
  # Final output
  p_summary <- list(
    p_Model_Array = p_Model_Array,
    p_RaosQ = p_RaosQ,
    p_RaosQ_rand = p_RaosQ_rand,
    NA_contrast = NA_contrast
  )
  
  return(p_summary)
}