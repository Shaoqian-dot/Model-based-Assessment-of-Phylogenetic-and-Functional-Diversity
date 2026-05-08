Power_Calculation <- function(p_values) {
  Power <- mean(p_values < 0.05, na.rm = TRUE)
  return(Power)
}