get_lrt_p <- function(fit1, fit0){
  
  if(inherits(fit1, "try-error") ||
     inherits(fit0, "try-error")) return(NA)
  
  out <- try(anova(fit0, fit1), silent = TRUE)
  
  if(inherits(out, "try-error")) return(NA)
  
  out$`Pr(>Chisq)`[2]
}

get_lrt_p_safe <- function(fit1, fit0){
  
  if (inherits(fit1, "try-error"))
    return(NA_real_)
  
  if (inherits(fit0, "try-error"))
    return(NA_real_)
  
  get_lrt_p(fit1, fit0)
}

