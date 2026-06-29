# When globalTest == TRUE, only poisson and gaussian families are allowed;
# When globalTest == FALSE and Swap == TRUE, only poisson and binomial families are allow; 
# When globalTest == FALSE and Swap == FALSE, only poisson and gaussian families are allow.
# run_simulation <- function(signal_type,
#                            family_type,
#                            globalTest,
#                            test,
#                            Swap, 
#                            p,
#                            q, 
#                            r,
#                            c_val,
#                            sigma2,
#                            seed,
#                            tree,
#                            alpha, 
#                            beta, 
#                            quantile, 
#                            Corr,
#                            Eigen,
#                            Method,
#                            Methods_para){ 
#   set.seed(seed)
#   
#   cat("Running:", signal_type,
#       "| family =", family_type,
#       "| p =", p, "\n")
#   
#   DM_phy_func <- getPhyloMatrix(tree = tree, m = p)
#   VC_phy_func <- 1 - DM_phy_func / max(DM_phy_func)
#   
#   I_p <- diag(p) 
#   J <- I_p - 1/p * matrix(1, p, p)   # centering matrix
#   S_J <- t(J) %*% VC_phy_func %*% J  # centered similarity matrix
#   eig <- spectral_decomp(VC_phy_func = S_J)
#   V_J <- eig$P # eigenvectors of the new similarity matrix
#   
#   if (Swap == TRUE) {
#     # Number of swap operations
#     NOPS <- switch(as.character(p),
#                    "5"  = 1,
#                    "10" = 2,
#                    "20" = 4,
#                    "40" = 8,
#                    "80" = 16,
#                    "160" = 32,
#                    stop("Unknown p")) 
#     dat <- getData(DM_phy_func = DM_phy_func, 
#                    NOPS = NOPS, 
#                    q = q, 
#                    alpha = alpha, 
#                    beta = beta, 
#                    r = r, 
#                    quantile = quantile,
#                    Corr = Corr, 
#                    P = V_J[, c(1 : (p - 1))], 
#                    Eigen = Eigen, 
#                    Distribution = family_type, 
#                    p = p) # Family_type can be "binomial" and "poisson" in swap simulations
#   dat_long <- dat$yX
#   dat_wide <- dat$abundance
#     } else {
#     dat <- generate_data(
#       p = p,
#       r = r,
#       signal = signal_type,
#       family = family_type,
#       c_val = c_val,
#       sigma2 = sigma2,
#       tree = tree,
#       V = V_J
#     )$long
#     dat_long <- dat$long
#     dat_wide <- dat$wide[, -1] # Remove the com column (the first column)
#   }
#   
#   model_res <- fit_models(
#     dat = dat_long,
#     family = family_type,
#     globalTest = globalTest,
#     tree = tree,
#     S_J = S_J,
#     p = p,
#     q = q,
#     test = test,
#     Method = Method,
#     Methods_para = Methods_para
#   ) 
#   
#   model_res <- model_res %>%
#     pivot_longer(
#       cols = -c(com, test),
#       names_to = "model",
#       values_to = "p_values"
#     ) %>%
#     mutate(method = "model_based")
#   
#   rao <- getRaosQ(
#     abundance = dat_wide,
#     DM_phy_func = DM_phy_func,
#     use_randomization = 0,
#     q = q
#   )
#   
#   rao_rand <- getRaosQ(
#     abundance = dat_wide,
#     DM_phy_func = DM_phy_func,
#     use_randomization = 1,
#     q = q
#   )
#   
#   rao_res <- bind_rows(
#     tibble(
#       com = 2:q,
#       model = "RaosQ",
#       p_values = rao
#     ),
#     tibble(
#       com = 2:q,
#       model = "Randomized RaosQ",
#       p_values = rao_rand
#     )
#   ) %>%
#     mutate(method = "raoQ")
#   
#   tmp_res <- bind_rows(model_res, 
#                        rao_res)
#   
#   res <- tmp_res %>%
#     mutate(
#       signal = signal_type,
#       family = family_type,
#       p = p,
#       r = r,
#       c_val = c_val,
#       globalTest = globalTest,
#       seed = seed,
#       Eigen = Eigen,
#       .before = 1
#     )
#   
#   return(res)
# }

run_simulation <- function(
    signal_type,
    family_type,
    globalTest,
    test,
    Swap,
    p,
    q,
    r,
    c_val,
    sigma2,
    seed,
    nsim,
    tree,
    alpha,
    beta,
    quantile,
    Corr,
    Eigen,
    Method,
    Methods_para
){
  res <- purrr::map_dfr(seq_len(nsim), function(sim){
    
    seed_i <- seed + sim - 1
    set.seed(seed_i)
    
    DM_phy_func <- getPhyloMatrix(tree = tree, m = p)
    VC_phy_func <- 1 - DM_phy_func / max(DM_phy_func)
    
    I_p <- diag(p)
    J <- I_p - 1 / p * matrix(1, p, p)
    S_J <- t(J) %*% VC_phy_func %*% J
    
    eig <- spectral_decomp(VC_phy_func = S_J)
    V_J <- eig$P
    
    cat(
      "Running:",
      signal_type,
      "| family =", family_type,
      "| p =", p,
      "| sim =", sim,
      "\n"
    )
    
    if (Swap) {
      
      NOPS <- switch(
        as.character(p),
        "5" = 1,
        "10" = 2,
        "20" = 4,
        "40" = 8,
        "80" = 16,
        "160" = 32,
        stop("Unknown p")
      )
      
      dat <- getData(
        DM_phy_func = DM_phy_func,
        NOPS = NOPS,
        q = q,
        alpha = alpha,
        beta = beta,
        r = r,
        quantile = quantile,
        Corr = Corr,
        P = V_J[, 1:(p - 1)],
        Eigen = Eigen,
        Distribution = family_type,
        p = p
      )
      
      dat_long <- dat$yX
      dat_wide <- dat$abundance
      
    } else {
      
      dat <- generate_data(
        p = p,
        r = r,
        signal = signal_type,
        family = family_type,
        c_val = c_val,
        sigma2 = sigma2,
        tree = tree,
        V = V_J
      )
      
      dat_long <- dat$long
      dat_wide <- dat$wide[, -1]
      
    }
    
    model_res <- fit_models(
      dat = dat_long,
      family = family_type,
      globalTest = globalTest,
      tree = tree,
      S_J = S_J,
      p = p,
      q = q,
      test = test,
      Method = Method,
      Methods_para = Methods_para
    )
    
    model_res <- model_res %>%
      tidyr::pivot_longer(
        cols = -c(com, test),
        names_to = "model",
        values_to = "p_values"
      ) %>%
      dplyr::mutate(method = "model_based")
    
    rao <- getRaosQ(
      abundance = dat_wide,
      DM_phy_func = DM_phy_func,
      use_randomization = 0,
      q = q
    )
    
    rao_rand <- getRaosQ(
      abundance = dat_wide,
      DM_phy_func = DM_phy_func,
      use_randomization = 1,
      q = q
    )
    
    rao_res <- dplyr::bind_rows(
      tibble::tibble(
        com = 2:q,
        model = "RaosQ",
        p_values = rao
      ),
      tibble::tibble(
        com = 2:q,
        model = "Randomized RaosQ",
        p_values = rao_rand
      )
    ) %>%
      dplyr::mutate(method = "raoQ")
    
    tmp_res <- dplyr::bind_rows(
      model_res,
      rao_res
    )
    
    tmp_res %>%
      dplyr::mutate(
        seed = seed_i,
        signal = signal_type,
        family = family_type,
        p = p,
        r = r,
        c_val = c_val,
        globalTest = globalTest,
        Eigen = Eigen,
        .before = 1
      )
    
  })
  
  return(res)
}