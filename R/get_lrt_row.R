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
  
  lrt_no_rr <- get_lrt_p_safe(
    fit_glmm_no_rr$fit,
    fit_glmm_no_rr_null$fit
  )
  
  lrt_rr <- get_lrt_p_safe(
    fit_glmm_rr$fit,
    fit_glmm_rr_null$fit
  )
  
  data.frame(
    com = j,
    glmm_no_rr = lrt_no_rr$p_value,
    glmm_no_rr_status = lrt_no_rr$status,
    glmm_no_rr_warning_msgs = fit_glmm_no_rr$warning_msgs,
    glmm_no_rr_null_warning_msgs = fit_glmm_no_rr_null$warning_msgs,
    glmm_no_rr_convergence = fit_glmm_no_rr$convergence,
    glmm_no_rr_pdHess = fit_glmm_no_rr$pdHess,
    glmm_rr = lrt_rr$p_value,
    glmm_rr_status = lrt_rr$status,
    glmm_rr_warning_msgs = fit_glmm_rr$warning_msgs,
    glmm_rr_null_warning_msgs = fit_glmm_rr_null$warning_msgs,
    glmm_rr_convergence = fit_glmm_rr$convergence,
    glmm_rr_pdHess = fit_glmm_rr$pdHess
  )
  
  # data.frame(
  #   com = j,
  #   glmm_no_rr = get_lrt_p_safe(
  #     fit_glmm_no_rr,
  #     fit_glmm_no_rr_null
  #   ),
  #   glmm_rr = get_lrt_p_safe(
  #     fit_glmm_rr,
  #     fit_glmm_rr_null
  #   )
  # )
}