# If globalTest == TRUE, only "LRT" test is applicable. 
fit_models <- function(dat, 
                       family = c("gaussian", "poisson", "binomial"), 
                       globalTest, 
                       tree, 
                       test = c("LRT", "Wald", "Both"),
                       S_J,
                       p,
                       q,
                       axis_method = c('Method1', 'Method2', 'Both')){
  # Likelihood Ratio Test
  get_lrt_row <- function(j, fit_glmm_no_rr, fit_glmm_rr, val_num.eig){
    
    com_left <- setdiff(2:q, j)
    
    if(length(com_left) == 0)
      com_left <- NULL
    
    fit_glmm_no_rr_null <- try(
      fit_model(
        type = "p_mid",
        matrix_type = "P_J",
        rr = FALSE,
        yX = dat_long,
        P_J = V_J,
        val_num.eig = val_num.eig,
        Distribution = fam,
        null = TRUE,
        com_left = com_left,
        p = p
      ),
      silent = TRUE
    )
    
    fit_glmm_rr_null <- try(
      fit_model(
        type = "p_mid",
        matrix_type = "P_J",
        rr = TRUE,
        yX = dat_long,
        P_J = V_J,
        val_num.eig = val_num.eig,
        Distribution = fam,
        null = TRUE,
        com_left = com_left,
        p = p
      ),
      silent = TRUE
    )
    
    data.frame(
      com = j,
      glmm_no_rr = get_lrt_p_safe(
        fit_glmm_no_rr,
        fit_glmm_no_rr_null
      ),
      glmm_rr = get_lrt_p_safe(
        fit_glmm_rr,
        fit_glmm_rr_null
      )
    )
  }
  
  # Wald Test
  get_wald_row <- function(j, fit_glmm_no_rr, fit_glmm_rr, val_num.eig){
    
    label <- paste0(
      "p",
      seq_len(val_num.eig),
      ":com",
      j
    )
    
    chi_square_no_rr <- getChisquare_fixed(
      Model = fit_glmm_no_rr,
      label = label
    )
    
    chi_square_rr <- getChisquare_fixed(
      Model = fit_glmm_rr,
      label = label
    )
    
    data.frame(
      com = j,
      glmm_no_rr = pchisq(
        chi_square_no_rr,
        df = val_num.eig,
        lower.tail = FALSE
      ),
      glmm_rr = pchisq(
        chi_square_rr,
        df = val_num.eig,
        lower.tail = FALSE
      )
    )
  }
  
  family <- match.arg(family)
  #test <- match.arg(test)
  
  if(test == "Both"){
    test <- c("Wald", "LRT")
  } else {
    
  }
  test <- match.arg(
    test,
    choices = c("Wald", "LRT"),
    several.ok = TRUE
  )
  
  dat_long <- dat
  
  fam <- switch(
    family,
    gaussian = gaussian(),
    poisson  = nbinom2,
    binomial = binomial()
  )
  
  if (globalTest == TRUE){
    fit_glmm_no_rr <- try(
      glmmTMB(
        y ~ com * sp,
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
    
    fit_glmm_no_rr_null <- try(
      glmmTMB(
        y ~ sp + com,
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
    
    fit_glmm_rr <- try(
      glmmTMB(
        y ~ com * sp +
          rr(sp + 0 | id, d = 2),
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
    
    fit_glmm_rr_null <- try(
      glmmTMB(
        y ~ sp + com +
          rr(sp + 0 | id, d = 2),
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
    return(data.frame(
      com = "global",
      glmm_no_rr = get_lrt_p_safe(fit_glmm_no_rr, fit_glmm_no_rr_null),
      glmm_rr    = get_lrt_p_safe(fit_glmm_rr, fit_glmm_rr_null)
    ))
  } else if (globalTest == FALSE) {
    
    axis_method_tbl <- data.frame(
      axis_method = c("Method1", "Method3"),
      Methods_para = c(2, NA)
    )
    
    if (axis_method == "Both") {
      axis_method_tbl <- axis_method_tbl
    } else {
      axis_method_tbl <- subset(axis_method_tbl, axis_method == axis_method)
    }
    
    eig <- spectral_decomp(VC_phy_func = S_J)
    V_J <- eig$P # eigenvectors of the new similarity matrix
    D_J <- eig$D # eigenvalues of the new similarity matrix
    
    # res_all <- pmap(axis_method_tbl, function(axis_method, Methods_para) {
    #   
    #   val_num.eig <- get_num.eig(
    #     Methods = axis_method,
    #     Methods_para = Methods_para,
    #     D = D_J
    #   )
    #   
    #   fit_glmm_no_rr <- try(
    #     fit_model(
    #       type = "p_mid",
    #       matrix_type = "P_J",
    #       rr = FALSE,
    #       yX = dat_long,
    #       P_J = V_J,
    #       val_num.eig = val_num.eig,
    #       Distribution = fam,
    #       null = FALSE,
    #       p = p
    #     ),
    #     silent = TRUE
    #   )
    #   
    #   fit_glmm_rr <- try(
    #     fit_model(
    #       type = "p_mid",
    #       matrix_type = "P_J",
    #       rr = TRUE,
    #       yX = dat_long,
    #       P_J = V_J,
    #       val_num.eig = val_num.eig,
    #       Distribution = fam,
    #       null = FALSE,
    #       p = p
    #     ),
    #     silent = TRUE
    #   )
    #   
    #   res_list <- list()
    #   
    #   if ("LRT" %in% test) {
    #     res_list$LRT <- do.call(
    #       rbind,
    #       lapply(2:q, 
    #              get_lrt_row,
    #              val_num.eig = val_num.eig,
    #              fit_glmm_no_rr = fit_glmm_no_rr,
    #              fit_glmm_rr = fit_glmm_rr)
    #     )
    #   }
    #   
    #   if ("Wald" %in% test) {
    #     res_list$Wald <- do.call(
    #       rbind,
    #       lapply(2:q, 
    #              get_wald_row,
    #              val_num.eig = val_num.eig,
    #              fit_glmm_no_rr = fit_glmm_no_rr,
    #              fit_glmm_rr = fit_glmm_rr)
    #     )
    #   }
    #   
    #   bind_rows(
    #     lapply(names(res_list), function(x) {
    #       cbind(
    #         axis_method = axis_method,
    #         test = x,
    #         res_list[[x]]
    #       )
    #     })
    #   )
    #   
    # })
    
    res <- purrr::pmap_dfr(
      grid,
      function(p, r, family_type, Eigen, sim, seed_i) {
        
        set.seed(seed_i)
        
        dir.create("timing_logs", showWarnings = FALSE)
        timing_file <- file.path(
          "timing_logs",
          sprintf("timing_seed_%d.csv", seed)
        )
        
        timing <- system.time({
          
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
            "| r =", r,
            "| Eigen =", Eigen,
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
            axis_method = axis_method
          )
          
          model_res <- model_res %>%
            tidyr::pivot_longer(
              cols = -c(axis_method, com, test),
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
          
          result <- dplyr::bind_rows(
            model_res,
            rao_res
          ) %>%
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
          
        })  # end system.time()
        
        timing_df <- tibble::tibble(
          seed = seed_i,
          signal = signal_type,
          family = family_type,
          p = p,
          r = r,
          Eigen = Eigen,
          elapsed_sec = timing["elapsed"]
        )
        
        write.table(
          timing_df,
          file = timing_file,
          sep = ",",
          row.names = FALSE,
          col.names = !file.exists(timing_file),
          append = TRUE
        )
        
        cat(
          sprintf(
            "Finished: family=%s p=%d r=%d Eigen=%d | %.2f hours\n",
            family_type,
            p,
            r,
            Eigen,
            timing["elapsed"] / 3600
          )
        )
        
        result
      }
    )
    
    res <- bind_rows(res_all)
    
    # val_num.eig <- get_num.eig(Methods = Method, Methods_para = Methods_para, D = D_J)
    # fit_glmm_no_rr <- try(
    #   fit_model(type = 'p_mid', matrix_type = 'P_J', rr = FALSE,
    #                              yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
    #                              Distribution = fam, null = FALSE, p = p),
    #   silent = TRUE
    # )
    # 
    # fit_glmm_rr <- try(
    #   fit_model(type = 'p_mid', matrix_type = 'P_J', rr = TRUE,
    #             yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
    #             Distribution = fam, null = FALSE, p = p),
    #   silent = TRUE
    # )
    # 
    # res_list <- list()
    # 
    # if ("LRT" %in% test) {
    #   
    #   res_list$LRT <- do.call(
    #     rbind,
    #     lapply(2:q, get_lrt_row)
    #   )
    #   
    # }
    # 
    # if ("Wald" %in% test) {
    #   
    #   res_list$Wald <- do.call(
    #     rbind,
    #     lapply(2:q, get_wald_row)
    #   )
    #   
    # }
    # res <- bind_rows(
    #   lapply(
    #     names(res_list),
    #     function(x) {
    #       cbind(
    #         test = x,
    #         res_list[[x]]
    #       )
    #     }
    #   )
    # )
    
    return(res)
    }
  }

