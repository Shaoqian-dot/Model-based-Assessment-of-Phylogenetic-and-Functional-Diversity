getChisquare_fixed <- function(Model, label){
  
  if (inherits(Model, "glmmTMB")){
    
   # Extract coefficient estimates (conditional model)
    coe_Model <-  summary(Model)$coefficients$cond
    
    # Extract variance-covariance matrix (conditional model)
    vcov_Model <- vcov(Model)$cond
    
    # Check if label exists
    if (!all(label %in% rownames(coe_Model))) {
      stop("Label not found in model coefficients")
    }
    
    # Retrieve the estimate and variance for the specified label
    est <- coe_Model[label, 'Estimate']
    var <- vcov_Model[label, label]
    
    # Compute chi-square statistic using the provided function
    chi_square <- getChisquare(var_L_beta = var,  coe_L_beta = est)
  } else if (length(Model) == 1 && is.na(Model)) {
    
    chi_square <- NA # Model is NA because of singular contrast matrix 
    
  } else {
    stop('Non-identifiable model object')
  }
  
  return(chi_square)
} 