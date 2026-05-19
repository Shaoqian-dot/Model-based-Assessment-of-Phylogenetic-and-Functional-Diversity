##################################################
## Main simulation runner
##################################################

families <- c("gaussian", "poisson")
signals <- c("v1", "vp")

res_all <- expand.grid(family = families,
                       signal = signals,
                       stringsAsFactors = FALSE) |>
  pmap_dfr(function(family, signal){
    
    run_simulation(signal, family)
  })