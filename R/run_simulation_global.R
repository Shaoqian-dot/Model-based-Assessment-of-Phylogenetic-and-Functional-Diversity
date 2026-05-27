run_simulation <- function(signal_type,
                           family_type,
                           globalTest,
                           p,
                           r,
                           c_val,
                           sigma2,
                           nsim,
                           seed,
                           tree){
  set.seed(seed)
  
  cat("Running:", signal_type,
      "| family =", family_type,
      "| p =", p, "\n")
  
  tmp_res <- map_dfr(
    1:nsim,
    function(sim){
      
      dat <- generate_data(
        p = p,
        r = r,
        signal = signal_type,
        family = family_type,
        c_val = c_val,
        sigma2 = sigma2,
        tree = tree
      )
      
      fit_models(dat, family = family_type, globalTest, tree = tree)
    }
  )
  
  power_res <- colMeans(tmp_res < 0.05, na.rm = TRUE)
  
  res <- data.frame(
    signal = signal_type,
    family = family_type,
    p = p,
    c_val = c_val,
    model = names(power_res),
    power = as.numeric(power_res),
    globalTest = globalTest,
    seed = seed
  )
  
  return(res)
}