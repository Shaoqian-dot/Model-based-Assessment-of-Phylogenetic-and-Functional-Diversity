getChisquare_random <- function (Model, val_num.eig, P, com) {
  rf <- ranef(Model, condVar = TRUE)
  
  # Posterior variances of the random slope (community effect) across species
  VC_beta_hat_ori <- attr(rf$cond$sp, "condVar")[com, com, ]
  
  # Conditional modes of random slopes
  Rcoef_beta_hat <- ranef(Model)$cond$sp[, com]
  
  # Apply sum-to-zero contrast to random effects
  #Rcoef_beta_hat_C <- OT(Rcoef_beta_hat, diag(VC_beta_hat_ori), length(Rcoef_beta_hat))$beta_hat_C
  #VC_beta_hat_C    <- OT(Rcoef_beta_hat, diag(VC_beta_hat_ori), length(Rcoef_beta_hat))$VC_beta_hat_C
  VC_beta_hat <- diag(VC_beta_hat_ori)
  coe_L_beta <-
    t(P[, 1:val_num.eig]) %*% Rcoef_beta_hat
  
  var_L_beta <-
    t(P[, 1:val_num.eig]) %*% VC_beta_hat %*% P[, 1:val_num.eig]
  
  chi_square <- getChisquare(var_L_beta, coe_L_beta)
  return(chi_square)
} 