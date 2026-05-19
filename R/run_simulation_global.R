run_simulation <- function(signal_type, family_type){
  
  res_all <- list()
  
  counter <- 1
  
  for(p in p_vec){
    
    cat("Running:", signal_type, "| family =", family_type, "| p =", p, "\n")
    
    tmp_res <- map_dfr(
      1:nsim,
      function(sim){
        
        dat <- generate_data(
          p = p,
          r = r,
          signal = signal_type,
          family = family_type,
          c_val = c_val,
          sigma2 = sigma2
        )
        
        fit_models(dat, family = family_type)
      }
    )
    
    power_res <- colMeans(tmp_res < 0.05, na.rm = TRUE)
    
    res_all[[counter]] <- data.frame(
      signal = signal_type,
      family = family_type,
      p = p,
      model = names(power_res),
      power = as.numeric(power_res)
    )
    
    counter <- counter + 1
  }
  
  bind_rows(res_all)
}
