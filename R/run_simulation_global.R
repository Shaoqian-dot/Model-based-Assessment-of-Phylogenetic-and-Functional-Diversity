# When globalTest == TRUE, only poisson and gaussian families are allowed;
# When globalTest == FALSE and Swap == TRUE, only poisson and binomial families are allow; 
# When globalTest == FALSE and Swap == FALSE, only poisson and gaussian families are allow.
run_simulation <- function(signal_type,
                           family_type,
                           globalTest,
                           test,
                           Swap, 
                           p,
                           q, 
                           r,
                           c_val,
                           sigma2,
                           nsim,
                           seed,
                           tree,
                           alpha, 
                           beta, 
                           quantile, 
                           Corr,
                           Eigen){ 
  set.seed(seed)
  
  cat("Running:", signal_type,
      "| family =", family_type,
      "| p =", p, "\n")
  
  tmp_res <- map_dfr(
    1:nsim,
    function(sim){
      DM_phy_func <- getPhyloMatrix(tree = tree, m = p)
      VC_phy_func <- 1 - DM_phy_func / (max(DM_phy_func) + 1)
      V <- spectral_decomp(VC_phy_func)$P 
      if (Swap == TRUE) {
        # Number of swap operations
        NOPS <- switch(as.character(p),
                       "5"  = 1,
                       "10" = 2,
                       "20" = 4,
                       "40" = 8,
                       "80" = 16,
                       "160" = 32,
                       stop("Unknown p")) 
        dat <- getData(DM_phy_func = DM_phy_func, NOPS = NOPS, q = q, 
                       alpha = alpha, beta = beta, r = r, quantile = quantile,
                       Corr = Corr, P = V, Eigen = Eigen, Distribution = family_type)$yX # Family_type can be "binomial" and "poisson" in swap simulations
      } else {
        dat <- generate_data(
          p = p,
          r = r,
          signal = signal_type,
          family = family_type,
          c_val = c_val,
          sigma2 = sigma2,
          tree = tree,
          V = V
        )$long
      }
      fit_models(dat, 
                 family = family_type, 
                 globalTest, 
                 tree = tree, 
                 Phy_SM = VC_phy_func, 
                 p = p, 
                 q = q,
                 test = test)
    }
  )
  
  power_res <- tmp_res %>%
    group_by(com) %>%
    summarise(
      across(
        where(is.numeric),
        ~ mean(.x < 0.05, na.rm = TRUE)
      ),
      .groups = "drop"
    )
  
  res <- power_res %>%
    pivot_longer(
      cols = -com,
      names_to = "model",
      values_to = "power"
    ) %>%
    mutate(
      signal = signal_type,
      family = family_type,
      p = p,
      c_val = c_val,
      globalTest = globalTest,
      seed = seed,
      .before = 1
    )
  
  # if (globalTest == TRUE) {
  #   power_res <- colMeans(tmp_res < 0.05, na.rm = TRUE)
  #   
  #   res <- data.frame(
  #     signal = signal_type,
  #     family = family_type,
  #     p = p,
  #     c_val = c_val,
  #     model = names(power_res),
  #     power = as.numeric(power_res),
  #     globalTest = globalTest,
  #     seed = seed
  #   )
  # } else {
  #   power_res <- tmp_res %>%
  #     group_by(com) %>%
  #     summarise(
  #       across(
  #         -com,
  #         ~ mean(.x < 0.05, na.rm = TRUE)
  #       ),
  #       .groups = "drop"
  #     )
  #   
  #   res <- power_res %>%
  #     pivot_longer(
  #       cols = -com,
  #       names_to = "model",
  #       values_to = "power"
  #     ) %>%
  #     mutate(
  #       signal = signal_type,
  #       family = family_type,
  #       p = p,
  #       c_val = c_val,
  #       globalTest = globalTest,
  #       seed = seed,
  #       .before = 1
  #     )
  # }
  return(res)
}