fit_model <- function(type, matrix_type, rr,
                      yX, P, P_J, val_num.eig, Distribution, 
                      null, com_left, p, site = FALSE, LV_num = 2){
  # Function for model fitting
  fit_once <- function(optCtrl = NULL, 
                       method = NULL, jitter.sd = 0) {
    msgs <- character()
    
    fit <- withCallingHandlers(
      glmmTMB(
        form,
        family = Distribution,
        data = yX_tmp,
        REML = FALSE,
        control = glmmTMBControl(
          start_method = list(method = method, jitter.sd = jitter.sd),
          optCtrl = optCtrl
        )
      ),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
      }
    )
    
    list(fit = fit, warnings = unique(msgs))
  }
  
  # Function for dealing with warning containing iteration limit or function evaluation limit reached
  fit_with_limit_retry <- function(attempt_seed, method = NULL, jitter.sd = 0) {
    set.seed(attempt_seed)
    result <- fit_once(method = method, jitter.sd = jitter.sd)
    
    limit_warning <- any(grepl(
      "iteration limit reached without convergence|function evaluation limit reached without convergence",
      result$warnings,
      ignore.case = TRUE
    ))
    
    if (limit_warning) {
      # Reuse the same randomized start; change only the optimizer limits.
      set.seed(attempt_seed)
      result <- tryCatch(
        fit_once(list(iter.max = 1000, eval.max = 1000), 
                 method = method, jitter.sd = jitter.sd),
        error = function(e) {
          message("Higher-limit retry failed: ", conditionMessage(e))
          result
        }
      )
    }
    
    result
  }
  # acceptable() decides if a model is converged and pdHess == TRUE. 
  #If yes, acceptable() returns TRUE; if not, acceptable() returns FALSE.  
  acceptable <- function(result) {
    isTRUE(result$fit$fit$convergence == 0L) &&
      isTRUE(result$fit$sdr$pdHess) &&
      isTRUE(is.finite(result$fit$fit$objective))
  }
  
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
      method <- 'res'
      jitter.sd <- 0.1
      
      base_formula <- paste(
        base_formula,
        "+ rr(sp + 0 | id, ",LV_num,")"
      )
      form <- as.formula(base_formula)
      
      # Distinct seeds give the initial fit and each refit different randomized starts.
      attempt_seeds <- sample.int(.Machine$integer.max, 21L)
      
      selected <- fit_with_limit_retry(attempt_seeds[1L], 
                                       method = method, jitter.sd = jitter.sd) # Dealing with warnings contain limit iteration/evaluations
      selected_attempt <- 0L  # 0 means the initial fit
      n_accepted_refits <- 0L
      
      needs_refit <- !acceptable(selected)
      
      if (needs_refit & !is.null(method)) {
        best_refit <- NULL
        
        fits <- vector("list", 21)
        fits[[1]] <- selected
        for (i in 2 : 21) {
          candidate <- tryCatch(
            fit_with_limit_retry(attempt_seeds[i], 
                                 method = method, jitter.sd = jitter.sd),
            error = function(e) {
              message("Refit ", i, " failed: ", conditionMessage(e))
              NULL
            }
          )
          
          fits[[i]] <- candidate
          
          if (is.null(candidate) || !acceptable(candidate)) next # next indicates jumping into the next for - loop
          
          n_accepted_refits <- n_accepted_refits + 1L
          
          if (is.null(best_refit) ||
              candidate$fit$fit$objective < best_refit$fit$fit$objective) {
            best_refit <- candidate
            selected_attempt <- i
          }
        }
        
        # Select a best refit.
        negativeLogLik_group <- c()
        for (i in 1 :21) {
          negativeLogLik_group[i] <- fits[[i]]$fit$fit$objective
        }
        
        # The following is my strategy that select a converged model having relatively smaller negative log-likelihood;
        #otherwise we use the model having the smallest negative log likelihood no matter is converge or not.
        if(is.null(best_refit)) {
          selected <- fits[[which(negativeLogLik_group == min(negativeLogLik_group))]]
        } else if (abs(best_refit$fit$fit$objective - min(negativeLogLik_group)) < 1) {
          selected <- best_refit
        } else {
          selected <- fits[[which(negativeLogLik_group == min(negativeLogLik_group))]]
        }
        
        # if (!is.null(best_refit) & 
        #     abs(best_refit$fit$fit$objective - min(negativeLogLik_group)) < 1) {
        #   # 1 is the threshold we can tolerate between the selected model and the model having the smallest negative log-likelihood 
        #   selected <- best_refit
        # } else {
        #   selected <- fits[[which(negativeLogLik_group == min(negativeLogLik_group))]]
        # }
        #assign("fits", fits, envir = .GlobalEnv)
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
    } else {
      form <- as.formula(base_formula)
      selected <- fit_with_limit_retry(attempt_seed = NULL, 
                                       method = NULL, jitter.sd = 0) # Dealing with warnings contain limit iteration/evaluations
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
        selected_attempt = NA, # the time of refit, which is as the best model
        needs_refit = NA # A logical variable to indicate if we need refit or not.
      ))
    }
  }
}