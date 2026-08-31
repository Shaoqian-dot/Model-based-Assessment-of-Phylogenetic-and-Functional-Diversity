# -----------------------------------------------------------------------------
# Boxplot with Dunn's test and compact letter display
# -----------------------------------------------------------------------------
#
# Arguments:
#   data         : data frame containing group and value columns
#   group_levels : desired order of groups
#   ylab         : y-axis label
#
# Returns:
#   A ggplot object.

Boxplot <- function(
    data,
    group_levels = c(
      "LF",
      "MF",
      "SF",
      "NE",
      "SE",
      "CG",
      "OP",
      "OA"
    ),
    ylab = "Randomized Rao's Q"
) {
  
  # Set the order of groups.
  data$group <- factor(
    data$group,
    levels = group_levels
  )
  
  
  # ---------------------------------------------------------------------------
  # Dunn's test with Benjamini-Hochberg correction
  # ---------------------------------------------------------------------------
  
  stat.test <- data %>%
    rstatix::dunn_test(
      value ~ group,
      p.adjust.method = "BH"
    )
  
  
  # ---------------------------------------------------------------------------
  # Compact letter display
  # ---------------------------------------------------------------------------
  
  cld <- multcompView::multcompLetters(
    setNames(
      stat.test$p.adj,
      paste(
        stat.test$group1,
        stat.test$group2,
        sep = "-"
      )
    )
  )
  
  letters_df <- data.frame(
    group = names(cld$Letters),
    label = cld$Letters
  )
  
  # Ensure the same group order is used in the plot.
  letters_df$group <- factor(
    letters_df$group,
    levels = group_levels
  )
  
  
  # ---------------------------------------------------------------------------
  # Position significance letters above the boxes
  # ---------------------------------------------------------------------------
  
  y_pos <- max(
    data$value,
    na.rm = TRUE
  ) * 1.1
  
  
  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  
  fig <- ggplot(
    data,
    aes(
      x = group,
      y = value,
      fill = group
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      alpha = 0.7
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.5,
      size = 2
    ) +
    geom_text(
      data = letters_df,
      aes(
        x = group,
        y = y_pos,
        label = label
      ),
      inherit.aes = FALSE,
      size = 5
    ) +
    theme_classic() +
    theme(
      legend.position = "none"
    ) +
    labs(
      x = "Group",
      y = ylab
    )
  
  return(fig)
}