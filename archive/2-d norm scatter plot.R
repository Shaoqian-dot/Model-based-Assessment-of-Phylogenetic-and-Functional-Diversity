# 安装并加载需要的包
# install.packages(c("MASS", "ggplot2", "patchwork"))

library(MASS)
library(ggplot2)
library(patchwork)

set.seed(123)

# -----------------------------
# 生成两个二维正态分布样本
# -----------------------------

n <- 500

mu1 <- c(1, 1)
mu2 <- c(1.5, 1.5)

Sigma <- matrix(c(1, 0.8,
                  0.8, 1), nrow = 2)

# 样本
x1 <- MASS::mvrnorm(n, mu = mu1, Sigma = Sigma)
x2 <- MASS::mvrnorm(n, mu = mu2, Sigma = Sigma)

df1 <- data.frame(x = x1[,1],
                  y = x1[,2],
                  group = "Blue")

df2 <- data.frame(x = x2[,1],
                  y = x2[,2],
                  group = "Red")

df <- rbind(df1, df2)

# -----------------------------
# 中间散点图
# -----------------------------

p_scatter <- ggplot(df, aes(x = x, y = y, color = group)) +
  geom_point(alpha = 0.5, size = 1.8) +
  scale_color_manual(values = c("Blue" = "blue",
                                "Red" = "red")) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  ) +
  labs(x = "X", y = "Y")

# -----------------------------
# 上侧边际密度（X）
# -----------------------------

p_top <- ggplot(df, aes(x = x, fill = group, color = group)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  scale_fill_manual(values = c("Blue" = "blue",
                               "Red" = "red")) +
  scale_color_manual(values = c("Blue" = "blue",
                                "Red" = "red")) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 0, 5.5)
  ) +
  labs(y = "Density")

# -----------------------------
# 右侧边际密度（Y）
# -----------------------------

p_right <- ggplot(df, aes(x = y, fill = group, color = group)) +
  geom_density(alpha = 0.25, linewidth = 1) +
  coord_flip() +
  scale_fill_manual(values = c("Blue" = "blue",
                               "Red" = "red")) +
  scale_color_manual(values = c("Blue" = "blue",
                                "Red" = "red")) +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 5.5, 0)
  ) +
  labs(x = "Density")

# -----------------------------
# 拼图
# -----------------------------

layout <- "
AB
CD
"

final_plot <-
  p_top + plot_spacer() +
  p_scatter + p_right +
  plot_layout(
    design = layout,
    widths = c(4, 1.2),
    heights = c(1.2, 4)
  )

final_plot
