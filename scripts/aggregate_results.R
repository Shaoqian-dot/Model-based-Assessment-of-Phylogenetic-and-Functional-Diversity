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
  here("output", "Diversity_Poisson_p20_Both_Method1_2_S2_RE(p-1-k)"),
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
    r,
    c_val,
    model,
    globalTest,
    test,
    Eigen,
    com,
    method
  ) %>%
  summarise(
    power = mean(p_values < 0.05, na.rm = TRUE),
    n_jobs = n(),
    .groups = "drop"
  )

############################################################
### Power plot
############################################################

# ======================
# User settings
# ======================

x_var <- "r"      # 改成 "p" 即可画 Power vs p
facet_var <- "com"

# ======================
# Data filter
# ======================

plot_data <- res_summary %>%
  filter(
    family == "poisson",
    Eigen == 2,
    globalTest == FALSE,
    com %in% c(2, 3, 4, 5)
  ) %>%
  mutate(
    test_plot = ifelse(is.na(test), "Baseline", test)
  )

x_breaks <- sort(unique(plot_data[[x_var]]))

# ======================
# Plot
# ======================

fig <- ggplot(
  plot_data,
  aes(
    x = .data[[x_var]],
    y = power,
    color = model,
    linetype = test_plot,
    group = interaction(model, test_plot)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(
    ~ com,
    nrow = 1
  ) +
  scale_linetype_manual(
    values = c(
      Wald = "dashed",
      LRT = "solid",
      Baseline = "dotdash"
    )
  ) +
  scale_color_manual(
    values = c(
      glmm_no_rr = "#1B9E77",
      glmm_rr = "#D95F02",
      RaosQ = "#7570B3",
      `Randomized RaosQ` = "#E7298A"
    )
  ) +
  scale_x_continuous(
    breaks = x_breaks
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  theme_bw() +
  labs(
    x = ifelse(
      x_var == "r",
      "Sample size (r)",
      "Number of species (p)"
    ),
    y = "Power",
    color = "Model",
    linetype = "Test"
  )
fig

############################################################
## Create results directory
############################################################
path <- 'figures/p20/poisson/Method1'
dir.create(
  here(path),
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
## Save figures
############################################################
ggsave(
  filename = here(path, "Diversity_Poisson_p20_Eigen2.png"),
  plot = fig,
  width = 16,
  height = 4.5,
  dpi = 300
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