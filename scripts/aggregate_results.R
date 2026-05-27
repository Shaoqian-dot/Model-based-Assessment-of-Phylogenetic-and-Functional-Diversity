rm(list = ls())

############################################################
## Aggregate simulation results
############################################################

library(dplyr)
library(purrr)
library(here)

cat("Aggregating simulation outputs...\n")

############################################################
## Locate output files
############################################################

files <- list.files(
  here("output"),
  pattern = "\\.rds$",
  full.names = TRUE
)

cat("Number of files found:", length(files), "\n")

if(length(files) == 0){
  stop("No .rds files found in output/")
}

############################################################
## Read all result files
############################################################

res_all <- map_dfr(files, readRDS)

############################################################
## Aggregate power estimates
############################################################

res_summary <- res_all %>%
  group_by(
    family,
    signal,
    p,
    c_val,
    model
  ) %>%
  summarise(
    power = mean(power, na.rm = TRUE),
    n_jobs = n(),
    .groups = "drop"
  )

############################################################
## Create results directory
############################################################

dir.create(
  here("results"),
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
## Save aggregated results
############################################################

saveRDS(
  res_summary,
  here("results", "power_summary.rds")
)

write.csv(
  res_summary,
  here("results", "power_summary.csv"),
  row.names = FALSE
)

############################################################
## Finish
############################################################

cat("Aggregation completed.\n")
cat("Results saved to results/\n")