fit_model <- function(type, matrix_type, rr,
                      yX, P, P_J, val_num.eig, Distribution){
  
  # Number of species
  m <- ncol(P)
  q <- 5
  # Select eigenvector matrix
  P_use <<- if (matrix_type == "P") P else P_J###############################################################
  
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
  
  Lsp <- Sp %*% P_use[, 1:m]
  colnames(Lsp) <- paste0("p", 1:m)
  
  yX_tmp <- cbind(yX, Lsp)
  
  p_vars_1 <- paste0("p", 1:val_num.eig)
  p_vars_2 <- paste0("p", (val_num.eig + 1) : m)
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
  colnames(Lsp) <- paste0("p", 1:m)

  yX_tmp <<- cbind(yX, Lsp)

  form_1 <- as.formula(
    paste("~ 0 + com * (",fixed_part_2,")")
  )
  X_propto <- model.matrix(form_1, data = yX_tmp)
  # n_random <- ncol(X_propto) # the number of random effects in 0 + com * fixed_part_2
  # K <- diag( n_random)
  # colnames(K) <- colnames(X_propto)
  # colnames(K) <- colnames(X_propto)
  K <- diag(m - val_num.eig)
  colnames(K) <- paste0("p", (val_num.eig + 1) : m)
  rownames(K) <- paste0("p", (val_num.eig + 1) : m)
  
  if (type == "p_mid") {
    # base_formula <- paste(
    #   "y ~ com * (", fixed_part_1, ") + propto(0 + com * (", fixed_part_2, ") |rep(1, nrow(yX_tmp)), K)"
    # )
    if (null) {
      base_formula <- paste(
        "y ~ com + (", fixed_part_1, ") + propto(0 + ", fixed_part_2, " | com, K)"
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
      data = yX_tmp
    )
    return(fit)
  }
}