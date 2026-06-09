rm(list = ls())

############################################################
## Aggregate simulation results
############################################################

library(dplyr)
library(purrr)
library(here)
library(ggplot2)

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
    model,
    globalTest
  ) %>%
  summarise(
    power = mean(power, na.rm = TRUE),
    n_jobs = n(),
    .groups = "drop"
  )

############################################################
### Power plot
############################################################
res_summary_filter <- res_summary %>%
  filter(
    family == "poisson",
    c_val == 0.9,
    globalTest == FALSE
  )
p_vec <- unique(res_summary_filter$p)
ggplot(
  res_summary_filter,
  aes(x = p,
      y = power,
      color = model)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ signal) +
  scale_x_continuous(
    breaks = p_vec
  ) +
  ylim(0, 1) +
  theme_bw() +
  labs(
    x = "Number of species (p)",
    y = "Power",
    color = "Model"
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