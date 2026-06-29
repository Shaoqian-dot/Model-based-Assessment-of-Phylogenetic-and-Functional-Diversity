fit_model <- function(type, matrix_type, rr,
                      yX, P, P_J, val_num.eig, Distribution, null, com_left, p){
  # Select eigenvector matrix
  P_use <- if (matrix_type == "P") P else P_J###############################################################
  
  # ===== CASE 1: contrast-based model =====
  if (type == "contrast") {
    
    C <- as.matrix(P_use[, 1:val_num.eig])
    yX_tmp <- yX
    
    res <- tryCatch({
      contrasts(yX_tmp$sp) <- C
      yX_tmp <- nameContrast(yX_tmp, val_num.eig)
      
      if (!rr) {
        fit <- glmmTMB(
          y ~ 0 + sp + com + sp:com,
          family = Distribution,
          data = yX_tmp
        )
      } else {
        fit <- glmmTMB(
          y ~ 0 + sp + com + sp:com +
            rr(sp + 0 | id, 2),
          family = Distribution,
          data = yX_tmp
        )
      }
      fit
    }, error = function(e){
      message("contrast model error: ", e$message)
      NA
    })
    
    return(res)
  }
  
  # ===== CASE 2: eigenvector (p) model =====
  Sp <- model.matrix(~0 + sp, data = yX)
  # Lsp <- Sp %*% P_use[, 1:val_num.eig]
  # colnames(Lsp) <- paste0("p", 1:val_num.eig)
  
  Lsp <- Sp %*% P_use[, 1:p]
  colnames(Lsp) <- paste0("p", 1:p)
  
  yX_tmp <- cbind(yX, Lsp)
  
  p_vars_1 <- paste0("p", 1:val_num.eig)
  p_vars_2 <- paste0("p", (val_num.eig + 1) : p)
  fixed_part_1 <- paste(p_vars_1, collapse = " + ")
  fixed_part_2 <- paste(p_vars_2, collapse = " + ")
  
  if (type == "p") {
     base_formula <- paste(
       "y ~ com * (", fixed_part_1, ") + diag(0 + com |sp)"
     ) 
    # base_formula <- paste(
    #   "y ~ com * (", fixed_part_1, ") + diag(0 + com * (", fixed_part_2 ,") |sp)"
    # )

    if (rr) {
      base_formula <- paste(
        base_formula,
        "+ rr(sp + 0 | id, 2)"
      )
    }

    form <- as.formula(base_formula)

    fit <- glmmTMB(
      form,
      family = Distribution,
      data = yX_tmp
    )
    return(fit)
  }
  
  # ===== CASE 3: eigenvector (p_fix) model =====
  if (type == "p_fix") {

    base_formula <- paste(
      "y ~ com * (", fixed_part_1, ")"
    )

    if (rr) {
      base_formula <- paste(
        base_formula,
        "+ rr(sp + 0 | id, 2)"
      )
    }

    form <- as.formula(base_formula)

    fit <- glmmTMB(
      form,
      family = Distribution,
      data = yX_tmp
    )
    return(fit)
  }
#   
  # ===== CASE 4: eigenvector (p_mid) model =====
  Sp <- model.matrix(~0 + sp, data = yX)
  Lsp <- Sp %*% P_use
  colnames(Lsp) <- paste0("p", 1:p)

  yX_tmp <- cbind(yX, Lsp)

  K <- diag(p - val_num.eig)
  colnames(K) <- paste0("p", (val_num.eig + 1) : p)
  # rownames(K) <- paste0("p", (val_num.eig + 1) : (p - 1))
  
  if (type == "p_mid") {
    # base_formula <- paste(
    #   "y ~ com * (", fixed_part_1, ") + propto(0 + com * (", fixed_part_2, ") |rep(1, nrow(yX_tmp)), K)"
    # )
    if (null) {
      for (l in seq_along(com_left)) {
        yX_tmp[[paste0("com", com_left[l])]] <- 
          as.numeric(yX_tmp$com == com_left[l])
      }
      
      interaction_terms <- if(length(com_left) == 0){
        NULL
      } else {
        paste0("com", com_left, ":(", fixed_part_1, ")")
      }
      
      rhs_terms <- c(
        "com",
        paste0("(", fixed_part_1, ")"),
        interaction_terms,
        paste0("propto(0 + ", fixed_part_2, " | com, K)")
      )
      
      base_formula <- paste(
        "y ~",
        paste(rhs_terms, collapse = " + ")
      )
      
      } else {
      base_formula <- paste(
        "y ~ com * (", fixed_part_1, ") + propto(0 + ", fixed_part_2, " | com, K)"
      )
    }
    # base_formula <- paste(
    #   "y ~ com * (", fixed_part_1, ") + diag(", fixed_part_2, "|com)"
    # )
    if (rr) {
      base_formula <- paste(
        base_formula,
        "+ rr(sp + 0 | id, 2)"
      )
    }

    form <- as.formula(base_formula)

    fit <- glmmTMB(
      form,
      family = Distribution,
      data = yX_tmp,
      REML=FALSE
    )
    return(fit)
  }
}