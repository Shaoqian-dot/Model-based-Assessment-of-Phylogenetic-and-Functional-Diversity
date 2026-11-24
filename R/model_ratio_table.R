model_ratio_table <- function(p, r, res) {
  d <- res[which(
    res$p == p & res$r == r & res$model == "glmm_rr"
  ), , drop = FALSE]
  
  pick <- function(test, com) {
    d[which(d$test == test & d$com == com), , drop = FALSE]
  }
  
  alt  <- pick("LRT", 2)
  wald <- pick("Wald", 2)
  n <- length(unique(alt$seed))
  
  #if (n == 0L) stop("No matching data sets for this p and r.")
  
  # check_rows <- function(x) {
  #   if (nrow(x) != n || !setequal(x$seed, alt$seed)) {
  #     stop("Expected one row per data set; check for missing or duplicate rows.")
  #   }
  # }
  # 
  # check_rows(alt)
  # check_rows(wald)
  
  # NA counts as not meeting the condition; denominator remains all n data sets.
  ratio <- function(x, target) sum(x %in% target) / n
  
  out <- data.frame(
    `Model descriptions` = c(
      "Convergence = 0",
      "pdHess = TRUE",
      "Status = 'Success' (Wald)",
      "Status = 'Success' (LRT)",
      "Refit = FALSE"
    ),
    Alt = c(
      ratio(alt$alt_convergence, 0),
      ratio(alt$alt_pdHess, TRUE),
      ratio(wald$status, "Success"),
      NA_real_,
      ratio(alt$alt_needs_refit, FALSE)
    ),
    check.names = FALSE
  )
  
  for (com in 2:5) {
    null <- pick("LRT", com)
    #check_rows(null)
    
    out[[paste0("Null_com", com)]] <- c(
      ratio(null$null_convergence, 0),
      ratio(null$null_pdHess, TRUE),
      NA_real_,
      ratio(null$status, "Success"),
      ratio(null$null_needs_refit, FALSE)
    )
  }
  
  attr(out, "n_datasets") <- n
  out
}