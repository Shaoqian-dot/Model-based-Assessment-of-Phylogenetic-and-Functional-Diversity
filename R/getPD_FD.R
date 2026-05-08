getPD_FD <- function(Model, yX, VC_phy_func, val_num.eig, Contrast, m, P, K){
  var_coe_L_beta <- getVar_Coe_L_beta(Model = Model, K = K, Contrast = Contrast, 
                                      m = m, P = P, val_num.eig = val_num.eig, yX = yX)
  coe_L_beta <- var_coe_L_beta$coe_L_beta
  var_L_beta <- var_coe_L_beta$var_L_beta
  chi_square <- getChisquare(var_L_beta = var_L_beta, coe_L_beta = coe_L_beta)# Too many parameters in this function. The first three are enough. 
  p_value <- 1 - pchisq(chi_square, val_num.eig)
  return(list(p_value = p_value,
              chi_square = chi_square))
} # p_value is a value, which can change with different contrasts.