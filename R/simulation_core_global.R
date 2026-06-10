fit_models <- function(dat, family = c("gaussian", "poisson", "binomial"), 
                       globalTest, tree, test = c("LRT", "Wald"),
                       Phy_SM){
  
  family <- match.arg(family)
  test <- match.arg(test)
  
  dat_long <- dat$long
  
  fam <- if(family == "gaussian"){
    gaussian()
  } else if (family == "poisson") {
    nbinom2
  } else if (family == "binomial") {
    binomial()
  }
  
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
  } else {
    p <- nlevels(dat$long$sp)
    q <- nlevels(dat$long$com)
    I_p <- diag(p) 
    J <- I_p - 1/p * matrix(1, p, p)   # centering matrix
    S_J <- t(J) %*% Phy_SM %*% J  # centered similarity matrix
    V_J <- spectral_decomp(VC_phy_func = S_J)$P  # eigenvectors of transformed matrix
    D <- spectral_decomp(VC_phy_func = S_J)$D
    val_num.eig <- get_num.eig(Methods = 'Method3', Methods_para = NA, D = D)
    
    fit_glmm_no_rr <- try(
      fit_model(type = 'p_mid', matrix_type = 'P_J', rr = FALSE,
                                 yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
                                 Distribution = fam, null = FALSE),
      silent = TRUE
    )
    
    fit_glmm_rr <- try(
      fit_model(type = 'p_mid', matrix_type = 'P_J', rr = TRUE,
                yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
                Distribution = fam, null = FALSE),
      silent = TRUE
    )
    
    if (test == "LRT") {
      res <- do.call(
        rbind,
        lapply(2:q, function(j) {
          
          if (q == 2) {
            com_left <- NULL
          } else {
            com_left <- setdiff(2:q, j)
          }
          
          fit_glmm_no_rr_null <- try(
            fit_model(type = "p_mid", matrix_type = "P_J", rr = FALSE,
                      yX = dat_long, P_J = V_J,
                      val_num.eig = val_num.eig,
                      Distribution = fam,
                      null = TRUE, com_left = com_left),
            silent = TRUE
          )
          
          fit_glmm_rr_null <- try(
            fit_model(type = "p_mid", matrix_type = "P_J", rr = TRUE,
                      yX = dat_long, P_J = V_J,
                      val_num.eig = val_num.eig,
                      Distribution = fam,
                      null = TRUE, com_left = com_left),
            silent = TRUE
          )
          
          data.frame(
            com = j,
            glmm_no_rr = get_lrt_p(fit_glmm_no_rr, fit_glmm_no_rr_null),
            glmm_rr    = get_lrt_p(fit_glmm_rr, fit_glmm_rr_null)
          )
        })
      )
    } else if (test == "Wald") {
      res <- vector("list", q - 1)
      
      for (j in 2:q) {
        
        label <- paste0("com", j, ":p", 1:val_num.eig)
        
        chi_square_no_rr <- getChisquare_fixed(
          Model = fit_glmm_no_rr,
          label = label
        )
        
        chi_square_rr <- getChisquare_fixed(
          Model = fit_glmm_rr,
          label = label
        )
        
        res[[j - 1]] <- data.frame(
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
      
      res <- do.call(rbind, res)
    }
    
    # 
    # for(j in 2 : q){
    #   
    #   if (q == 2) {
    #     com_left <- NULL
    #   } else {
    #     com_left <- setdiff(2:q, j)
    #   }
    #   
    #   fit_glmm_no_rr_null <- try(
    #     fit_model(type = 'p_mid', matrix_type = 'P_J', rr = FALSE,
    #               yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
    #               Distribution = fam, null = TRUE, com_left = com_left),
    #     silent = TRUE
    #   )
    #   
    #   fit_glmm_rr_null <- try(
    #     fit_model(type = 'p_mid', matrix_type = 'P_J', rr = TRUE,
    #               yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
    #               Distribution = fam, null = TRUE, com_left = com_left),
    #     silent = TRUE
    #   )
    #   
    #   data.frame(
    #     glmm_no_rr = get_lrt_p(fit_glmm_no_rr, fit_glmm_no_rr_null),
    #     glmm_rr    = get_lrt_p(fit_glmm_rr, fit_glmm_rr_null)
    #   )
      
      # Model_p_mid_J_null <- fit_model(type = 'p_mid', matrix_type = 'P_J', rr = FALSE,
      #                                 yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, 
      #                                 Distribution = Distribution, null = TRUE, com_left = com_left)
      # 
      # Model_p_mid_J_rr_null <- fit_model(type = 'p_mid', matrix_type = 'P_J', rr = TRUE,
      #                                    yX = yX, P = P, P_J = P_J, val_num.eig = val_num.eig, 
      #                                    Distribution = Distribution, null = TRUE, com_left = com_left)
      # Model_p_mid_J_alt <- Model_p_mid_J
      # Model_p_mid_J_rr_alt <- Model_p_mid_J_rr
      # 
      # p_LR_Model_p_mid_J    <- anova(Model_p_mid_J_null, Model_p_mid_J_alt)$`Pr(>Chisq)`[2]
      # p_LR_Model_p_mid_J_rr <- anova(Model_p_mid_J_rr_null, Model_p_mid_J_rr_alt)$`Pr(>Chisq)`[2]
      # p_Model_Matrix[j - 1, ] <- c(p_LR_Model_p_mid_J, p_LR_Model_p_mid_J_rr)
    }
  }

