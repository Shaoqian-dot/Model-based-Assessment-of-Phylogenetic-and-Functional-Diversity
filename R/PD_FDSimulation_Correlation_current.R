
PD_FDSimulation <- function (alpha, beta, q, m, m_less, NOPS, r, diversity,
                             Methods, Methods_para, tree, trait, quantile,
                             Model = Model, Corr = Corr, n_Model, K = 10){ # num.eig can be replaced with `method`.
  # browser()
  # Initialization
  p_summary <- list()
   p_Model_Array <- array(NA, dim = c(q - 1, n_Model, length(Methods)))
   dimnames(p_Model_Array) <- list(
     row    = paste0("com", 2:q),
     col    = paste0("Model", 1:n_Model),
     face   =  c(
       "p1, p2",
       "λ > mean"
     )
   )
   
   {
      if (diversity == 1) {
        DM_phy_func <- getPhyloMatrix(tree = tree, m = m)$DM_phy_func
              }
      else if (diversity == 2) {
        DM_phy_func <- getFuncMatrix(trait = trait, m = m, m_less = m_less)
      }
    }
    VC_phy_func <- max(DM_phy_func) - DM_phy_func
    P <- spectral_decomp(VC_phy_func = VC_phy_func)$P # Eigenvectors matrix
    D <- spectral_decomp(VC_phy_func = VC_phy_func)$D
    num.eig <- get_num.eig(Methods = Methods, Methods_para = Methods_para, D = D)
    Data <- getData(DM_phy_func = DM_phy_func, NOPS = NOPS,
                    m = m, q = q, alpha = alpha, 
                    beta = beta, r = r, quantile = quantile,
                    Model = Model, Corr = Corr)#Contain randomness
    yX <- Data$yX
    abundance <- Data$abundance
    
    for (k in 1 : length(num.eig)) {
      val_num.eig <- num.eig[k]
      ################################################################################
      ## ===== MODEL FITTING =========================================================
      ################################################################################
      # Model 1: Eigenvector subset used as species contrasts
      C <- as.matrix(P[, c(1 : val_num.eig)])
      #rownames(C) <- levels(yX$sp)
      #colnames(C) <- c(paste0("C", 1 : val_num.eig))
      yX_C4 <- yX
      res <- tryCatch(
        {
          contrasts(yX_C4$sp) <- C
          yX_C4 <- nameContrast(m, yX_C4, val_num.eig)
          fit_glmm_fixed_C4 <- glmmTMB(
            y ~ 0 + sp + com + sp:com,
            family = binomial,
            data   = yX_C4
          )
          # Model 2: Model 1 with additional reduced-rank structure
          fit_glmm_fixed_C4_rr <- glmmTMB(
            y ~ 0 + sp + com + sp:com +
              rr(sp + 0 | id , 2),
            family = binomial,
            data   = yX_C4
          )
          list(fit_glmm_fixed_C4 = fit_glmm_fixed_C4,
               fit_glmm_fixed_C4_rr = fit_glmm_fixed_C4_rr)
        },
        error = function(e){
          message("Model 1 2 - Error occurred: ", e$message)
          list(fit_glmm_fixed_C4 = NA,
               fit_glmm_fixed_C4_rr = NA)
        }
      )
      fit_glmm_fixed_C4 <- res$fit_glmm_fixed_C4
      fit_glmm_fixed_C4_rr <- res$fit_glmm_fixed_C4_rr
      
      # Prepare data for Models 3–4 (David’s eigenvector-based formulation)
      Sp = model.matrix(~0+sp, data = yX)
      Lsp = Sp %*% P[, c(1 : val_num.eig)]
      colnames(Lsp) <- paste0('p', 1 : val_num.eig)
      yX_p_added <- cbind(yX, Lsp)
      
      p_vars <- paste0("p", 1:val_num.eig)
      fixed_part <- paste(p_vars, collapse = " + ")
      
      # Model 3: Eigenvector coefficients with random effects
      form_random <- as.formula(
        paste("y ~  com * (", fixed_part, ") + diag(0 + com|sp)")
      )
      fit_glmm_David_random <- glmmTMB(
        form_random,
        family = binomial,
        data = yX_p_added
      )
      
      # Model 4: Same as Model 3 with additional reduced-rank term
      form_random_rr <- as.formula(
        paste("y ~  com * (", fixed_part, ") + diag(0 + com|sp) + rr(sp + 0 | id, 2)")
      )
      fit_glmm_David_random_rr <- glmmTMB(
        form_random_rr,
        family = binomial,
        data = yX_p_added
      )
      
      # Model 5: random effect model
      fit_glmm_rf <- glmmTMB(
        y ~ 0 + sp + com + diag(com | sp),
        family = binomial,
        data   = yX
      )
      
      # Model 6 with rr
      fit_glmm_rf_rr <- glmmTMB(
        y ~ 0 + sp + com + diag(com | sp) +
          rr(sp + 0 | id, 2),
        family = binomial,
        data   = yX
      )
      
      # Model 7: New Similarity matrix: P_new as contrasts
      I_m <- diag(m) 
      J <- I_m - 1/m * matrix(1, m, m)
      S_new <- t(J) %*% VC_phy_func %*% J
      P_new <- spectral_decomp(VC_phy_func = S_new)$P
      C5 <- as.matrix(P_new [, c(1 : val_num.eig)])
      #rownames(C5) <- levels(yX$sp)
      #colnames(C5) <- c(paste0("C", 1 : val_num.eig))
      yX_C5 <- yX
      res2 <- tryCatch(
        {
          contrasts(yX_C5$sp) <- C5
          yX_C5 <- nameContrast(m, yX_C5, val_num.eig)
          
          fit_glmm_C5 <- glmmTMB(
            y ~ 0 + sp + com + sp:com,
            family = binomial,
            data   = yX_C5
          )
          # Model 8: Model 7 with additional reduced-rank structure
          fit_glmm_C5_rr <- glmmTMB(
            y ~ 0 + sp + com + sp:com +
              rr(sp + 0 | id , 2),
            family = binomial,
            data   = yX_C5
          )
          list(fit_glmm_C5 = fit_glmm_C5,
               fit_glmm_C5_rr = fit_glmm_C5_rr)
        },
        error = function(e){
          message("Model 7 8 - Error occurred: ", e$message)
          list(fit_glmm_C5 = NA,
               fit_glmm_C5_rr = NA)
        }
      )
      fit_glmm_C5 <- res2$fit_glmm_C5
      fit_glmm_C5_rr <- res2$fit_glmm_C5_rr
      
      # Matrix to store chi-square statistics
      chi_square <- matrix(NA, q - 1, n_Model)
      dimnames(chi_square) <- list(
        row    = paste0("com", 2:q),
        col    = paste0("Model", 1:n_Model)
      )
      for(j in 2 : q){
        label_Model_1278 <- paste0('spC', (1 : val_num.eig), ':com', j)
        label_Model_34 <- paste0('com', j, ':p', 1 : val_num.eig) 
        
        # Model 1
        chi_square_1 <- getChisquare_fixed (fit_glmm_fixed_C4, label_Model_1278)
        
        # Model 2
        chi_square_2 <- getChisquare_fixed (fit_glmm_fixed_C4_rr, label_Model_1278)
        
        # Model 3
        
        chi_square_3 <- getChisquare_fixed (fit_glmm_David_random, label_Model_34)
        
        # Model 4
        chi_square_4 <- getChisquare_fixed (fit_glmm_David_random_rr, label_Model_34)
        
        # Model 5
        chi_square_5 <- getChisquare_random(fit_glmm_rf, val_num.eig, P_new, j)
        
        # Model 6
        chi_square_6 <- getChisquare_random(fit_glmm_rf_rr, val_num.eig, P_new, j)
        
        # Model 7
        chi_square_7 <- getChisquare_fixed (fit_glmm_C5, label_Model_1278)
        
        # Model 8
        chi_square_8 <- getChisquare_fixed (fit_glmm_C5_rr, label_Model_1278)
        
        # Store chi-square statistics
        chi_square[j - 1, ] <- c(
          chi_square_1, chi_square_2, 
          chi_square_3, chi_square_4, 
          chi_square_5, chi_square_6,
          chi_square_7, chi_square_8
        )
      }
      p_Model_Matrix <- 1 - pchisq(chi_square, df = val_num.eig)
      p_Model_Array [, , k] <- p_Model_Matrix
      
      rm(fit_glmm_fixed_C4, fit_glmm_fixed_C4_rr, 
         fit_glmm_David_random, fit_glmm_David_random_rr,
         fit_glmm_rf, fit_glmm_rf_rr,
         fit_glmm_C5, fit_glmm_C5_rr)
      gc()
    }
    p_RaosQ <- getRaosQ(abundance = abundance, DM_phy_func = DM_phy_func,
                        randomization = 0, q = q, r = r, m = m)
    p_RaosQ_rand <- getRaosQ(abundance = abundance, DM_phy_func = DM_phy_func,
                                      randomization = 1, q = q, r = r, m = m)
    p_summary <- list(
      p_Model_Array = p_Model_Array,
      p_RaosQ = p_RaosQ,
      p_RaosQ_rand = p_RaosQ_rand
    )
  return(p_summary)
  }