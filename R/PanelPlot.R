PanelPlot <- function (x_var = "r",
                       q_fixed = 5,
                       eigen_target = 3,
                       axis_method_label = "Method3",
                       levels_p = c(5, 10, 20, 40, 80),
                       levels_r = c(4, 8, 16, 32, 64),
                       Methods_para,
                       res_summary) {
  row_var <- if (x_var == "r") "p" else "r"
  scenario_levels <- c(
    "Small ΔPD\nwith ΔComposition",
    "Large ΔPD\nwith ΔComposition",
    "No ΔPD\nwith ΔComposition",
    "No ΔPD\nwith ΔOverall abundance"
  )
  
  # Values displayed on the x-axis
  x_breaks <- if (x_var == "r") {
    q_fixed * levels_r
  } else {
    levels_p
  }
  
  # plot_data <- res_summary %>%
  #   filter(
  #     .data$family == "binomial",
  #     .data$Eigen == .env$eigen_target,
  #     ifelse(.env$axis_method_label == 'Method3',
  #            is.na(.data$axis_method) |
  #              .data$axis_method == .env$axis_method_label,
  #            is.na(.data$axis_method) |
  #              .data$axis_method == .env$axis_method_label
  #            & near(.data$Methods_para, .env$Methods_para))
  #     
  #   )
  plot_data <- res_summary %>%
    dplyr::filter(
      .data$family == "binomial",
      .data$Eigen == .env$eigen_target,
      if (.env$axis_method_label == "Method3") {
        is.na(.data$axis_method) |
          .data$axis_method == .env$axis_method_label
      } else {
        is.na(.data$axis_method) |
          (
            .data$axis_method == .env$axis_method_label &
              dplyr::near(.data$Methods_para, .env$Methods_para)
          )
      }
    ) %>%
    mutate(
      test_plot = coalesce(as.character(.data$test), "Baseline"),
      
      scenario = case_when(
        .data$com == 3 ~ "Small ΔPD\nwith ΔComposition",
        .data$com == 4 ~ "Large ΔPD\nwith ΔComposition",
        .data$com == 2 ~ "No ΔPD\nwith ΔComposition",
        .data$com == 5 ~ "No ΔPD\nwith ΔOverall abundance",
        TRUE ~ NA_character_
      ),
      
      scenario = factor(
        scenario,
        levels = scenario_levels
      ),
      
      # Actual variable plotted on the x-axis
      x_value = if (.env$x_var == "r") {
        .env$q_fixed * .data$r
      } else {
        .data$p
      },
      
      row_label = if (.env$row_var == "p") {
        factor(
          .data$p,
          levels = .env$levels_p,
          labels = paste("# Species =", .env$levels_p)
        )
      } else {
        factor(
          .data$r,
          levels = .env$levels_r,
          labels = paste(
            "Sample size =",
            .env$q_fixed * .env$levels_r
          )
        )
      }
    ) %>%
    filter(!is.na(.data$scenario))
  
  fig <- ggplot(
    plot_data,
    aes(
      x = x_value,
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
    
    scale_linetype_manual(
      values = c(
        Wald = "dashed",
        LRT = "solid",
        Baseline = "dotdash"
      )
    ) +
    
    scale_colour_manual(
      values = c(
        glmm_no_rr = "#1B9E77",
        glmm_rr = "#D95F02",
        RaosQ = "#7570B3",
        `Randomized RaosQ` = "#E7298A"
      )
    ) +
    
    scale_x_continuous(
      trans = "log2",
      breaks = x_breaks,
      labels = x_breaks
    ) +
    
    coord_cartesian(ylim = c(0, 1)) +
    
    labs(
      x = if (
        x_var == "r"
      ) {
        "Sample size (log2 scale)"
      } else {
        "Number of species (log2 scale)"
      },
      y = NULL,
      colour = "Model",
      linetype = "Test"
    ) +
    
    theme_bw() +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold"
      )
    )
  return(fig)
}