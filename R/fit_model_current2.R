fit_model <- function(type, matrix_type, rr,
                      yX, P, P_J, val_num.eig, Distribution, 
                      null, com_left, p, site = FALSE, LV_num = 2){
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
  
  Lsp <- Sp %*% P_use
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
  K <- diag(p - val_num.eig)
  colnames(K) <- paste0("p", (val_num.eig + 1) : p)
  
  if (type == "p_mid") {
    if (null) {
      for (l in seq_along(com_left)) {
        yX_tmp[[paste0("com", com_left[l])]] <- as.numeric(yX_tmp$com == com_left[l]) #Create community indicators
      }
      
      interaction_terms <- if(length(com_left) == 0){
        NULL
      } else {
        paste0("com", com_left, ":(", fixed_part_1, ")")
      }
      
      rhs_terms <- c(
        paste0("(", fixed_part_1, ")"),
        interaction_terms,
        paste0("propto(0 + ", fixed_part_2, " | com, K)")
      )
      
      base_formula <- paste(
        "y ~",
        paste(rhs_terms, collapse = " + ")
      )
      
      } else {
        if (site == TRUE) {
          base_formula <- paste(
            "y ~ (", fixed_part_1, ") + (1|site) + com : (", fixed_part_1, ") + propto(0 + ", fixed_part_2, " | com, K)"
          )
        } else {
          base_formula <- paste(
            "y ~ (", fixed_part_1, ") + com : (", fixed_part_1, ") + propto(0 + ", fixed_part_2, " | com, K)"
          )
        }
     }
    if (rr) {
      base_formula <- paste(
        base_formula,
        "+ rr(sp + 0 | id, ",LV_num,")"
      )
      form <- as.formula(base_formula)
      
      fit_once <- function(optCtrl = NULL) {
        msgs <- character()
        
        fit <- withCallingHandlers(
          glmmTMB(
            form,
            family = Distribution,
            data = yX_tmp,
            REML = FALSE,
            control = glmmTMBControl(
              start_method = list(method = "res", jitter.sd = 0.1),
              optCtrl = optCtrl
            )
          ),
          warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
          }
        )
        
        list(fit = fit, warnings = unique(msgs))
      }
    } else {
      form <- as.formula(base_formula)
      
      fit_once <- function(optCtrl = NULL) {
        msgs <- character()
        
        fit <- withCallingHandlers(
          glmmTMB(
            form,
            family = Distribution,
            data = yX_tmp,
            REML = FALSE,
            control = glmmTMBControl(
              optCtrl = optCtrl
            )
          ),
          warning = function(w) {
            msgs <<- c(msgs, conditionMessage(w))
          }
        )
        
        list(fit = fit, warnings = unique(msgs))
      }
    }

    fit_with_limit_retry <- function(attempt_seed) {
      set.seed(attempt_seed)
      result <- fit_once()
      
      limit_warning <- any(grepl(
        "iteration limit reached without convergence|function evaluation limit reached without convergence",
        result$warnings,
        ignore.case = TRUE
      ))
      
      if (limit_warning) {
        # Reuse the same randomized start; change only the optimizer limits.
        set.seed(attempt_seed)
        result <- tryCatch(
          fit_once(list(iter.max = 1000, eval.max = 1000)),
          error = function(e) {
            message("Higher-limit retry failed: ", conditionMessage(e))
            result
          }
        )
      }
      
      result
    }
    
    acceptable <- function(result) {
      isTRUE(result$fit$fit$convergence == 0L) &&
        isTRUE(result$fit$sdr$pdHess) &&
        isTRUE(is.finite(result$fit$fit$objective))
    }
    
    # Distinct seeds give the initial fit and each refit different randomized starts.
    attempt_seeds <- sample.int(.Machine$integer.max, 21L)
    
    selected <- fit_with_limit_retry(attempt_seeds[1L])
    selected_attempt <- 0L  # 0 means the initial fit
    n_accepted_refits <- 0L
    
    needs_refit <- !acceptable(selected)
    if (needs_refit) {
      best_refit <- NULL
      
      fits <- vector("list", 21)
      fits[[1]] <- selected
      for (i in seq_len(20L)) {
        candidate <- tryCatch(
          fit_with_limit_retry(attempt_seeds[i + 1L]),
          error = function(e) {
            message("Refit ", i, " failed: ", conditionMessage(e))
            NULL
          }
        )
        
        fits[[i + 1]] <- candidate
        
        if (is.null(candidate) || !acceptable(candidate)) next
        
        n_accepted_refits <- n_accepted_refits + 1L
        
        if (is.null(best_refit) ||
            candidate$fit$fit$objective < best_refit$fit$fit$objective) {
          best_refit <- candidate
          selected_attempt <- i
        }
        }
      
      # Replace selected with best_refit.
      if (!is.null(best_refit)) selected <- best_refit
      assign("fits", fits, envir = .GlobalEnv)
    }
    
    
    fit <- selected$fit #Note that if none of the 20 refits is acceptable, retain the initial fit.
    warning_msgs <- if (length(selected$warnings)) {
      paste(selected$warnings, collapse = " | ")
    } else {
      "None"
    }
    
    return(list(
      fit = fit,
      warning_msgs = warning_msgs,
      convergence = fit$fit$convergence,
      pdHess = fit$sdr$pdHess,
      selected_attempt = selected_attempt, # the time of refit, which is as the best model
      needs_refit = needs_refit # A logical variable to indicate if we need refit or not.
    ))
  }
}