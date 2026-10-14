# Wald Test
get_wald_row <- function(j, fit_glmm_no_rr, fit_glmm_rr, val_num.eig){
  
  label <- paste0(
    "p",
    seq_len(val_num.eig),
    ":com",
    j
  )
  
  chi_square_no_rr <- getChisquare_fixed(
    Model = fit_glmm_no_rr$fit,
    label = label
  )
  
  chi_square_rr <- getChisquare_fixed(
    Model = fit_glmm_rr$fit,
    label = label
  )
  
  data.frame(
    com = j,
    glmm_no_rr = pchisq(
      chi_square_no_rr,
      df = val_num.eig,
      lower.tail = FALSE
    ),
    glmm_no_rr_status = 'Success',
    glmm_no_rr_warning_msgs = fit_glmm_no_rr$warning_msgs,
    glmm_no_rr_convergence = fit_glmm_no_rr$convergence,
    glmm_no_rr_pdHess = fit_glmm_no_rr$pdHess,
    glmm_rr = pchisq(
      chi_square_rr,
      df = val_num.eig,
      lower.tail = FALSE
    ),
    glmm_rr_status = 'Success',
    glmm_rr_warning_msgs = fit_glmm_rr$warning_msgs,
    glmm_rr_convergence = fit_glmm_rr$convergence,
    glmm_rr_pdHess = fit_glmm_rr$pdHess
  )
}