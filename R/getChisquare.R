# Function to compute the chi-square statistic for a linear combination of coefficients
# var_L_beta: variance-covariance matrix of the linear combination
# coe_L_beta: vector of coefficients for the linear combination
getChisquare <- function(var_L_beta, coe_L_beta){
  
  # Ensure the coefficients are a column vector
  coe_L_beta <- as.vector(coe_L_beta)
  
  # Compute the chi-square statistic using the formula:
  # chi^2 = t(coe_L_beta) %*% solve(var_L_beta) %*% coe_L_beta
  chi_square <- coe_L_beta %*% solve(var_L_beta) %*% coe_L_beta
  
  # Return the chi-square value
  return(chi_square)
}