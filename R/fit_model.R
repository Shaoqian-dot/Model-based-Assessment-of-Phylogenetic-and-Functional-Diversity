fit_model <- function(type, matrix_type, rr,
                      yX, P, P_J, val_num.eig, Distribution){
  
  # Number of species
  m <- ncol(P)
  
  # Select eigenvector matrix
  P_use <- if (matrix_type == "P") P else P_J
  
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
  Lsp <- Sp %*% P_use[, 1:val_num.eig]
  colnames(Lsp) <- paste0("p", 1:val_num.eig)
  
  yX_tmp <- cbind(yX, Lsp)
  
  p_vars <- paste0("p", 1:val_num.eig)
  fixed_part <- paste(p_vars, collapse = " + ")
  
  if (type == "p") {
    
        
    base_formula <- paste(
      "y ~ com * (", fixed_part, ") + diag(0 + com|sp)"
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
  
  # ===== CASE 3: eigenvector (p_fix) model =====
  if (type == "p_fix") {

    base_formula <- paste(
      "y ~ com * (", fixed_part, ")"
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
  
  # ===== CASE 4: eigenvector (p) model =====
  Sp <- model.matrix(~0 + sp, data = yX)
  Lsp <- Sp %*% P_use
  colnames(Lsp) <- paste0("p", 1:m)
  
  yX_tmp <- cbind(yX, Lsp)
  
  p_vars_1 <- paste0("p", 1:val_num.eig)
  p_vars_2 <- paste0("p", (val_num.eig + 1) : m)
  fixed_part_1 <- paste(p_vars_1, collapse = " + ")
  fixed_part_2 <- paste(p_vars_2, collapse = " + ")
  
  if (type == "p_mid") {
    base_formula <- paste(
      "y ~ com * (", fixed_part_1, ") + propto(0 + ", fixed_part_2, " |com)"
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
}