rm(list = ls())
library(here)
library(dplyr)
source('R/PanelPlot.R')
source('R/model_ratio_table.R')
#folder <- 'Eigen1_binomial_axis_methodMethod3_res'
folder <- 'Phy_fam_binom_test_Both_axis_Method3_Eig_1_Abs_TRUE_p_1020_r_8163264_muRef_MS_default_S3_sim1000'

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

# Power and Type I error
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

PanelPlot (x_var = "r",
           q_fixed = 5,
           eigen_target = 1,
           axis_method_label = "Method3",
           levels_p <- c(10, 20),
           levels_r <- c(8, 16, 32, 64),
           res_summary = res_summary)

# Warning Summary
res = res_all
table_10_8 <- model_ratio_table(p = 20, r = 64, res = res)
print(table_10_8, digits = 3, row.names = FALSE)
attr(table_10_8, "n_datasets")  # 672
