rm(list = ls())
library(here)

folder <- 'Phy_fam_binomial_test_Both_axis_method_Method3_Eigen_3_Absence_TRUE_p_10_r_8_muRef_MS_Default'

files <- list.files(
  here("output", folder),
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

res_all <- purrr::map_dfr(files, function(f) {
  x <- readRDS(f)
  
  # Rename old column to the new name
  if ("p_values" %in% names(x)) {
    x <- dplyr::rename(x, pvalue = p_values)
  }
  
  x
})

res_summary <- res_all %>%
  group_by(
    family,
    p,
    r,
    model,
    test,
    Eigen,
    com,
    method,
    axis_method) %>%
  summarise(
    power = mean (pvalue < 0.05, na.rm = TRUE),
    n_jobs = n(),
    .groups = 'drop'
  )

res <- res_all
Target_Variable <- "status"
group_vars <- c(
  "p", "r", "axis_method", "model",
  "Eigen", "test", "com",
  Target_Variable
)

warning_summary <- res %>%
  group_by(across(all_of(group_vars))) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(
    p, r, axis_method, model,
    Eigen, test, com, desc(n)
  )


if(!is.null(best_refit) & 
   abs(best_refit$fit$fit$objective - 5) < 1){}
