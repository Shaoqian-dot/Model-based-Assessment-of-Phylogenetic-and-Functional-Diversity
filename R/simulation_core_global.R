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
    
    res_all <- pmap(axis_method_tbl, function(axis_method, Methods_para) {
      
      val_num.eig <- get_num.eig(
        Methods = axis_method,
        Methods_para = Methods_para,
        D = D_J
      )
      
      fit_glmm_no_rr <- try(
        fit_model(
          type = "p_mid",
          matrix_type = "P_J",
          rr = FALSE,
          yX = dat_long,
          P_J = V_J,
          val_num.eig = val_num.eig,
          Distribution = fam,
          null = FALSE,
          p = p
        ),
        silent = TRUE
      )
      
      fit_glmm_rr <- try(
        fit_model(
          type = "p_mid",
          matrix_type = "P_J",
          rr = TRUE,
          yX = dat_long,
          P_J = V_J,
          val_num.eig = val_num.eig,
          Distribution = fam,
          null = FALSE,
          p = p
        ),
        silent = TRUE
      )
      
      res_list <- list()
      
      if ("LRT" %in% test) {
        res_list$LRT <- do.call(
          rbind,
          lapply(2:q, 
                 get_lrt_row,
                 val_num.eig = val_num.eig,
                 fit_glmm_no_rr = fit_glmm_no_rr,
                 fit_glmm_rr = fit_glmm_rr)
        )
      }
      
      if ("Wald" %in% test) {
        res_list$Wald <- do.call(
          rbind,
          lapply(2:q, 
                 get_wald_row,
                 val_num.eig = val_num.eig,
                 fit_glmm_no_rr = fit_glmm_no_rr,
                 fit_glmm_rr = fit_glmm_rr)
        )
      }
      
      bind_rows(
        lapply(names(res_list), function(x) {
          cbind(
            axis_method = axis_method,
            test = x,
            res_list[[x]]
          )
        })
      )
      
    })
    
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

