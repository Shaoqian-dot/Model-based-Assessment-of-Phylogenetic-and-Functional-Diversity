#OA: Overall abundance
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
  here("output", "PhyD_Pois_p510204080_r_TestBoth_MethodBoth_Eigen2"),
  pattern = "\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

cat("Number of files found:", length(files), "\n")

if(length(files) == 0){
  stop("No .rds files found in output/")
}

############################################################
## Read all result files
############################################################

res_all <- purrr::map_dfr(files, function(f) {
  
  x <- readRDS(f)
  
  # Rename old column to the new name
  if ("p_values" %in% names(x)) {
    x <- dplyr::rename(x, pvalue = p_values)
  }
  
  x
})
############################################################
## Aggregate power estimates
############################################################

# res_summary <- res_all %>%
#   group_by(
#     family,
#     signal,
#     p,
#     r,
#     c_val,
#     model,
#     globalTest,
#     test,
#     Eigen,
#     com,
#     method,
#     axis_method
#   ) %>%
#   summarise(
#     power = mean(pvalue < 0.05, na.rm = TRUE),
#     n_jobs = n(),
#     .groups = "drop"
#   )

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
    method,
    axis_method
  ) %>%
  summarise(
    power = mean(pvalue < 0.05, na.rm = TRUE),
    alt_warning_ratio = mean(alt_warning_msgs != "None", na.rm = TRUE),
    null_warning_ratio = mean(null_warning_msgs != "None", na.rm = TRUE),
    pvalue_failure_ratio = mean(status != "Success"),
    n_jobs = n(),
    .groups = "drop"
  )

# ############################################################
# ### Power plot
# ############################################################
# 
# # ======================
# # User settings
# # ======================
# 
# x_var <- "r"      # 改成 "p" 即可画 Power vs p
# facet_var <- "com"
# 
# # ======================
# # Data filter
# # ======================
# 
# plot_data <- res_summary %>%
#   filter(
#     p == 10,
#     family == "poisson",
#     Eigen == 1,
#     globalTest == FALSE,
#     com %in% c(2, 3, 4, 5),
#     axis_method == "Method1" | is.na(axis_method)
#   ) %>%
#   mutate(
#     test_plot = ifelse(is.na(test), "Baseline", test)
#   )
# 
# x_breaks <- sort(unique(plot_data[[x_var]]))
# 
# # ======================
# # Plot
# # ======================
# 
# fig <- ggplot(
#   plot_data,
#   aes(
#     x = .data[[x_var]],
#     y = power,
#     color = model,
#     linetype = test_plot,
#     group = interaction(model, test_plot)
#   )
# ) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   facet_wrap(
#     ~ com,
#     nrow = 1
#   ) +
#   scale_linetype_manual(
#     values = c(
#       Wald = "dashed",
#       LRT = "solid",
#       Baseline = "dotdash"
#     )
#   ) +
#   scale_color_manual(
#     values = c(
#       glmm_no_rr = "#1B9E77",
#       glmm_rr = "#D95F02",
#       RaosQ = "#7570B3",
#       `Randomized RaosQ` = "#E7298A"
#     )
#   ) +
#   scale_x_continuous(
#     breaks = x_breaks
#   ) +
#   coord_cartesian(
#     ylim = c(0, 1)
#   ) +
#   theme_bw() +
#   labs(
#     x = ifelse(
#       x_var == "r",
#       "Sample size (r)",
#       "Number of species (p)"
#     ),
#     y = "Power",
#     color = "Model",
#     linetype = "Test"
#   )
# fig
#install.packages("ggh4x")
library(ggh4x)

x_var = 'p'
row_var <- ifelse(x_var == "r", "p", "r")         # choose "p" or "r" for panel rows
if (x_var == 'r') x_var <- "total_n"           # "r" or "p"
q <- 5
levels_p = c(5, 10, 20, 40, 80)
levels_r = c(4, 8, 16, 32, 64)
axis_method_label <- 'Method1'
plot_data <- res_summary %>%
  filter(
    family == "poisson",
    Eigen == 2,
    globalTest == FALSE,
    axis_method == axis_method_label | is.na(axis_method)
  ) %>%
  mutate(
    test_plot = ifelse(is.na(test), "Baseline", test),
    
    ## First two columns = Power
    ## Last two columns = Type I error
    metric = ifelse(com %in% c(2, 5), "Type I error", "Power"),
    
    ## Column labels
    scenario = dplyr::case_when(
      com == 3 ~ "Small ΔPD\nwith ΔComposition",
      com == 4 ~ "Large ΔPD\nwith ΔComposition",
      com == 2 ~ "No ΔPD\nwith ΔComposition",
      com == 5 ~ "No ΔPD\nwith ΔOverall abundace"
    )
  )

plot_data <- plot_data %>%
  mutate(
    row_label = if (row_var == "p") {
      factor(
        p,
        levels = levels_p,
        labels = paste("# Species =", levels_p)
      )
    } else {
      factor(
        r,
        levels = levels_r,
        labels = paste("Sample size =", q * levels_r)
      )
    }
  )
plot_data <- plot_data %>%
  mutate(total_n = q * r)

plot_data$scenario <- factor(
  plot_data$scenario,
  levels = c(
    "Small ΔPD\nwith ΔComposition",
    "Large ΔPD\nwith ΔComposition",
    "No ΔPD\nwith ΔComposition",
    "No ΔPD\nwith ΔOverall abundace"
  )
)

fig <- ggplot(
  plot_data,
  aes(
    x = .data[[x_var]],
    y = power,
    colour = model,
    linetype = test_plot,
    group = interaction(model, test_plot)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  ggh4x::facet_grid2(
    rows = vars(row_label),
    cols = vars(scenario),
    switch = "y"
  ) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.right = element_text(angle = 0, face = "bold")
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
    trans = "log2",
    breaks = if (x_var == "total_n") {
      q * levels_r
    } else {
      levels_p
    },
    labels = waiver()
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(
    x = ifelse(
      x_var == "total_n",
      "Sample size (log2 scale)",
      "Number of species (log2 scale)"
    ),
    y = NULL,
    colour = "Model",
    linetype = "Test"
  ) 

fig

############################################################
## Create results directory
############################################################
path <- "figures/poisson/Eigen2"

if (!dir.exists(here(path))) {
  dir.create(
    here(path),
    recursive = TRUE
  )
}

############################################################
## Save figures
############################################################
ggsave(
  filename = here(path, "Phy_Poisson_Eigen2_axisMethod1_Power#sp.png"),
  plot = fig,
  width = 8,
  height = 8,
  dpi = 300
)

# Plot for comparing eigenvector choosing methods 
x_var <- 'r'
row_var <- ifelse(x_var == "r", "p", "r")         # choose "p" or "r" for panel rows
if (x_var == 'r') x_var <- "total_n"           # "r" or "p"

plot_data <- res_summary %>%
  filter(
    family == "poisson",
    Eigen == 2,
    globalTest == FALSE,
    model == "glmm_rr",
    axis_method %in% c("Method1", "Method3"),
    com %in% c(2, 5)
  ) %>%
  mutate(
    scenario = case_when(
      com == 2 ~ "No ΔPD\nwith ΔComposition",
      com == 5 ~ "No ΔPD\nwith ΔOverall abundace"
    )
  )

plot_data <- plot_data %>%
  mutate(
    row_label = if (row_var == "p") {
      factor(
        p,
        levels = levels_p,
        labels = paste("# Species =", levels_p)
      )
    } else {
      factor(
        r,
        levels = levels_r,
        labels = paste("Sample size =", q * levels_r)
      )
    }
  )

plot_data <- plot_data %>%
  mutate(total_n = q * r)

plot_data$scenario <- factor(
  plot_data$scenario,
  levels = c(
    "No ΔPD\nwith ΔComposition",
    "No ΔPD\nwith ΔOverall abundace"
  )
)

fig <- ggplot(
  plot_data,
  aes(
    x = .data[[x_var]],
    y = power,
    colour = axis_method,
    linetype = test,
    group = interaction(axis_method, test)
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  ggh4x::facet_grid2(
    rows = vars(row_label),
    cols = vars(scenario),
    switch = "y"
  ) +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.right = element_text(
      angle = 0,
      face = "bold"
    )
  ) +
  scale_colour_manual(
    values = c(
      Method1 = "#1B9E77",
      Method3 = "#D95F02"
    )
  ) +
  scale_linetype_manual(
    values = c(
      LRT = "solid",
      Wald = "dashed"
    )
  ) +
  scale_x_continuous(
    trans = "log2",
    breaks = if (x_var == "total_n") q * levels_r else levels_p,
    labels = waiver()
  ) +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  labs(
    x = ifelse(
      x_var == "total_n",
      "Sample size (log2 scale)",
      "Number of species (log2 scale)"
    ),
    y = "Type I error",
    colour = "Axis method",
    linetype = "Test"
  )
fig

############################################################
## Create results directory
############################################################
path <- "figures/poisson/Eigen2"

if (!dir.exists(here(path))) {
  dir.create(
    here(path),
    recursive = TRUE
  )
}

############################################################
## Save figures
############################################################
ggsave(
  filename = here(path, "Phy_Poisson_Eigen2_Power_r_axisMethodscompare.png"),
  plot = fig,
  width = 5.5,
  height = 8,
  dpi = 300
)

############################################################
## Summarize warning messages
############################################################

warning_summary <- res_all %>%
  group_by(p, r, axis_method, model, Eigen, test, alt_warning_msgs,) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(p, r, axis_method, model, Eigen, test, desc(n))
warning_summary

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