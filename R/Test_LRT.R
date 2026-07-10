# get_lrt_p <- function(fit1, fit0){
#   
#   if(inherits(fit1, "try-error") ||
#      inherits(fit0, "try-error")) return(NA)
#   
#   out <- try(anova(fit0, fit1), silent = TRUE)
#   
#   if(inherits(out, "try-error")) return(NA)
#   
#   out$`Pr(>Chisq)`[2]
# }

# get_lrt_p <- function(fit1, fit0){
#   
#   if (inherits(fit1, "try-error") ||
#       inherits(fit0, "try-error"))
#     return(NA)
#   
#   out <- tryCatch(
#     anova(fit0, fit1),
#     error = function(e) e
#   )
#   
#   if (inherits(out, "error")) {
#     message("anova() failed: ", conditionMessage(out))
#     return(NA)
#   }
#   
#   out$`Pr(>Chisq)`[2]
# }
# 
# 
# get_lrt_p_safe <- function(fit1, fit0){
#   
#   if (inherits(fit1, "try-error"))
#     return(NA_real_)
#   
#   if (inherits(fit0, "try-error"))
#     return(NA_real_)
#   
#   get_lrt_p(fit1, fit0)
# }

get_lrt_p_safe <- function(fit1, fit0) {

  # Alternative model failed
  if (inherits(fit1, "try-error")) {
    return(list(
      p_value = NA_real_,
      status = "Alternative model failed"
    ))
  }

  # Null model failed
  if (inherits(fit0, "try-error")) {
    return(list(
      p_value = NA_real_,
      status = "Null model failed"
    ))
  }

  out <- try(anova(fit0, fit1), silent = TRUE)
  
  if (inherits(out, "try-error")) {
    return(list(
      p_value = NA_real_,
      status = "anova() failed"
    ))
  }
  
  ## Check the log-likelihoods
  # if (is.na(out$logLik[1])) {
  #   return(list(
  #     p_value = NA_real_,
  #     status = "Null model logLik NA"
  #   ))
  # }
  # 
  # if (is.na(out$logLik[2])) {
  #   return(list(
  #     p_value = NA_real_,
  #     status = "Alternative model logLik NA"
  #   ))
  # }

  ## Check the log-likelihoods
  null_logLik_na <- is.na(out$logLik[1])
  alt_logLik_na  <- is.na(out$logLik[2])
  
  if (null_logLik_na && alt_logLik_na) {
    return(list(
      p_value = NA_real_,
      status = "Both models logLik NA"
    ))
  }
  
  if (null_logLik_na) {
    return(list(
      p_value = NA_real_,
      status = "Null model logLik NA"
    ))
  }
  
  if (alt_logLik_na) {
    return(list(
      p_value = NA_real_,
      status = "Alternative model logLik NA"
    ))
  }
  # Successful
  list(
    p_value = out$`Pr(>Chisq)`[2],
    status = "Success"
  )
}



