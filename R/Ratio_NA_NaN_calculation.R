Ratio_NA_NaN_calculation <- function(p_values) {
  Ratio <- mean(is.na(p_values))
  return(Ratio)
}