rm(list = ls())
library(here)
library(dplyr)
library(ggplot2)
source('R/PanelPlot.R')
source('R/model_ratio_table.R')
# folder <- 'Phy_bino_Wald_Met3_Eig_1_Abse_TRUE_p_40_r_8163264_muRef_MS_default_S3_sim1000'
# #folder <- 'Phy_fam_bino_test_Both_axis_Method3_Eig_1_Abse_TRUE_p_1020_r_8163264_Ref_MS_res_S3_sim1000'
# 
# files <- list.files(
#   here("output", folder),
#   pattern = "\\.rds$",
#   recursive = TRUE,
#   full.names = TRUE
# )

folders <- c(
  #"output/Phy_binomial_Wald_Method3_Eigen_1_Absence_TRUE_p_80_r_8163264_muRef_MS_default_S3_sim1000",
  "output/Phy_bino_Wald_Met3_Eig_1_Abse_TRUE_p_40_r_8163264_muRef_MS_default_S3_sim1000"
  #"output/Phy_fam_bino_testBoth_Method3_Eigen_1_AbsenceTRUE_p_80_r_816_muRef_MS_res_S3_sim1000",
  #"output/Phy_fam_bino_Both_Method3_Eigen_1_Absence_TRUE_numFit3_p_80_r_32_muRef_MS_res_S3_sim1000"
)

files <- unlist(lapply(
  folders,
  list.files,
  pattern = "\\.rds$",
  full.names = TRUE
))

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
  dplyr::filter((test != "LRT" | is.na(test)) 
                & !(res_all$seed %in% unique(res_all$seed[is.na(res_all$pvalue)]))) %>%
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
           levels_p <- c(40),
           levels_r <- c(8, 16, 32, 64),
           res_summary = res_summary)

# Warning Summary
res = res_all
table_10_8 <- model_ratio_table(p = 80, r = 64, res = res)
print(table_10_8, digits = 3, row.names = FALSE)
attr(table_10_8, "n_datasets")  # 672
