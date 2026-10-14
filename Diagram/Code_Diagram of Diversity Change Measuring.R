# ##Old codes
# #install.packages("metR")
# # Load required libraries
# library(ggplot2)
# library(reshape2)
# library(grid)
# library(metR)
# 
# # Define a 2x2 positive-definite correlation matrix
# S <- matrix(c(1, 0.8,
#               0.8, 1), nrow = 2)
# 
# # Eigen decomposition
# eig <- eigen(S)
# P <- eig$vectors
# lambda <- eig$values
# 
# # Define function f(beta) = beta^T S beta
# f_beta <- function(b1, b2) {
#   beta <- c(b1, b2)
#   return(t(beta) %*% S %*% beta)
# }
# 
# # Create grid
# x <- seq(-2, 2, length.out = 300)
# y <- seq(-2, 2, length.out = 300)
# z <- outer(x, y, Vectorize(function(b1, b2) f_beta(b1, b2)))
# 
# # Prepare data for ggplot
# df <- melt(z)
# df$Var1 <- x[df$Var1]
# df$Var2 <- y[df$Var2]
# 
# # Scale eigenvectors for arrows
# arrow_len <- 1.2
# p1 <- P[,1] * arrow_len
# p2 <- P[,2] * arrow_len
# 
# # Generate unit circle
# theta <- seq(0, 2*pi, length.out = 300)
# circle_df <- data.frame(x = cos(theta), y = sin(theta))
# 
# # Create plot
# ggplot(df, aes(x = Var1, y = Var2, z = value)) +
#   # Contour lines with labels
#   geom_contour(aes(z = value), color = "black", bins = 20) +
#   geom_text_contour(aes(z = value), stroke = 0.2, skip = 0, check_overlap = TRUE) +
#   
#   # Eigenvector arrows
#   geom_segment(aes(x = 0, y = 0, xend = p1[1], yend = p1[2]),
#                arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
#   geom_segment(aes(x = 0, y = 0, xend = p2[1], yend = p2[2]),
#                arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
#   
#   # Eigenvector labels
#   annotate("text", x = p1[1]*1.1, y = p1[2]*1.1, label = expression(bold(p)[1]), color = "red", parse = TRUE, size = 5) +
#   annotate("text", x = p2[1]*1.1, y = p2[2]*1.1, label = expression(bold(p)[2]), color = "red", parse = TRUE, size = 5) +
#   
#   # Unit circle (dashed)
#   geom_path(data = circle_df, aes(x = x, y = y), linetype = "dashed", color = "blue", inherit.aes = FALSE) +
#   
#   # Formula annotation
#   annotate("text", x = -1.8, y = 1.8,
#            label = expression(f(beta) == beta^T * S * beta),
#            parse = TRUE, size = 5) +
#   
#   # Axis labels
#   labs(x = expression(beta[1]), y = expression(beta[2])) +
#   coord_fixed() +
#   theme_minimal(base_size = 14)
# ##Old codes
# 
# 
# 
# 
# 
# ##New code
# #install.packages("metR")
# library(ggplot2)
# library(reshape2)
# library(grid)
# library(metR)
# 
# S <- matrix(c(1, 0.8,
#               0.8, 1), nrow = 2)
# 
# eig <- eigen(S)
# P <- eig$vectors
# lambda <- eig$values
# 
# f_beta <- function(b1, b2) {
#   beta <- c(b1, b2)
#   return(t(beta) %*% S %*% beta)
# }
# 
# x <- seq(-2, 2, length.out = 300)
# y <- seq(-2, 2, length.out = 300)
# z <- outer(x, y, Vectorize(function(b1, b2) f_beta(b1, b2)))
# 
# df <- melt(z)
# df$Var1 <- x[df$Var1]
# df$Var2 <- y[df$Var2]
# 
# arrow_len <- 1.2
# p1 <- P[,1] * arrow_len
# p2 <- P[,2] * arrow_len
# 
# theta <- seq(0, 2*pi, length.out = 300)
# circle_df <- data.frame(x = cos(theta), y = sin(theta))
# 
# ggplot(df, aes(x = Var1, y = Var2, z = value)) +
#   geom_contour(aes(z = value), color = "black", bins = 20) +
#   
#   # Remove contour labels
#   # geom_text_contour(aes(z = value), stroke = 0.2, skip = 0, check_overlap = TRUE) +
#   
#   # Eigenvectors
#   geom_segment(aes(x = 0, y = 0, xend = p1[1], yend = p1[2]),
#                arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
#   geom_segment(aes(x = 0, y = 0, xend = p2[1], yend = p2[2]),
#                arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
#   annotate("text", x = p1[1]*1.1, y = p1[2]*1.1, label = expression(bold(p)[1]), color = "red", parse = TRUE, size = 5) +
#   annotate("text", x = p2[1]*1.1, y = p2[2]*1.1, label = expression(bold(p)[2]), color = "red", parse = TRUE, size = 5) +
#   
#   # Unit circle
#   geom_path(data = circle_df, aes(x = x, y = y), linetype = "dashed", color = "blue", inherit.aes = FALSE) +
#   
#   # Formula annotation
#   annotate("text", x = -1.8, y = 1.8,
#            label = expression(f(beta) == beta^T * S * beta),
#            parse = TRUE, size = 5) +
#   
#   # # Explanation text
#   # annotate("text", x = -1.8, y = -1.5, hjust = 0,
#   #          label = "β: beta-diversity change vector\nf: Phylogenetic diversity change\nS: similarity matrix\np1, p2: eigenvectors of S",
#   #          size = 4, color = "black") +
#   
#   # Axis and title
#   labs(x = expression(beta[1]), y = expression(beta[2])) +
#   ggtitle(expression("Contour plot of " ~ f(beta) == beta^T * S * beta)) +
#   coord_fixed() +
#   theme_minimal(base_size = 14)
# ##New code
# 


##New code 1 for making the text right
# install.packages("cowplot")
library(ggplot2)
library(reshape2)
library(grid)
library(metR)
library(cowplot)

# Define matrix and eigen decomposition
S <- matrix(c(1, 0.8,
              0.8, 1), nrow = 2)
eig <- eigen(S)
P <- eig$vectors
lambda <- eig$values

# Define f(beta)
f_beta <- function(b1, b2) {
  beta <- c(b1, b2)
  return(t(beta) %*% S %*% beta)
}

# Generate grid
x <- seq(-2, 2, length.out = 300)
y <- seq(-2, 2, length.out = 300)
z <- outer(x, y, Vectorize(function(b1, b2) f_beta(b1, b2)))
df <- melt(z)
df$Var1 <- x[df$Var1]
df$Var2 <- y[df$Var2]

# Arrows (eigenvectors)
arrow_len <- 1.2
p1 <- P[,1] * arrow_len
p2 <- P[,2] * arrow_len

# Unit circle
theta <- seq(0, 2*pi, length.out = 300)
circle_df <- data.frame(x = cos(theta), y = sin(theta))

# Main plot without explanation text
main_plot <- ggplot(df, aes(x = Var1, y = Var2, z = value)) +
  geom_contour(aes(z = value), color = "black", bins = 20) +
  geom_segment(aes(x = 0, y = 0, xend = p1[1], yend = p1[2]),
               arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
  geom_segment(aes(x = 0, y = 0, xend = p2[1], yend = p2[2]),
               arrow = arrow(length = unit(0.3, "cm")), color = "red", size = 1) +
  annotate("text", x = p1[1]*1.1, y = p1[2]*1.1, label = expression(bold(p)[1]), color = "red", parse = TRUE, size = 5) +
  annotate("text", x = p2[1]*1.1, y = p2[2]*1.1, label = expression(bold(p)[2]), color = "red", parse = TRUE, size = 5) +
  geom_path(data = circle_df, aes(x = x, y = y), linetype = "dashed", color = "blue", inherit.aes = FALSE) +
  annotate("text", x = -1.8, y = 1.8,
           label = expression(f(beta) == beta^T * S * beta),
           parse = TRUE, size = 5) +
  labs(x = expression(beta[1]), y = expression(beta[2])) +
  ggtitle(expression("Contour plot of " ~ f(beta) == beta^T * S * beta)) +
  coord_fixed() +
  theme_minimal(base_size = 14)

# Text-only plot for explanation
explanation_text <- ggdraw() + 
  draw_label("β: beta-diversity change vector\nf: Phylogenetic diversity change\nS: similarity matrix\np1, p2: eigenvectors of S", 
             hjust = 0, x = 0, y = 0.8, size = 12)

# Combine plots side by side
combined <- plot_grid(main_plot, explanation_text, nrow = 1, rel_widths = c(3, 1))

# Print the combined plot
print(combined)
##New code 2 for making the text right










library(ggplot2)
library(grid)

# Similarity matrix and eigendecomposition
S <- matrix(c(1, 0.8,
              0.8, 1),
            nrow = 2, byrow = TRUE)

eig <- eigen(S)

# Orient eigenvectors consistently
v1_unit <- eig$vectors[, 1]

if (sum(v1_unit) < 0) {
  v1_unit <- -v1_unit
}

# Orthogonal second eigenvector
v2_unit <- c(-v1_unit[2], v1_unit[1])

# Length of displayed eigenvector arrows
arrow_len <- 1.25
v1 <- arrow_len * v1_unit
v2 <- arrow_len * v2_unit

# Generate quadratic-form surface
x <- seq(-2, 2, length.out = 300)
y <- seq(-2, 2, length.out = 300)

df <- expand.grid(
  beta1 = x,
  beta2 = y
)

df$value <- with(
  df,
  S[1, 1] * beta1^2 +
    2 * S[1, 2] * beta1 * beta2 +
    S[2, 2] * beta2^2
)

# Unit circle
circle_angle <- seq(0, 2 * pi, length.out = 400)

circle_df <- data.frame(
  x = cos(circle_angle),
  y = sin(circle_angle)
)

# ------------------------------------------------------------
# Define community-change vectors
#
# com 2 is nearly perpendicular to v1:
#     small projection onto v1
#
# com 3 is close to v1:
#     large projection onto v1
# ------------------------------------------------------------

theta2 <- 80 * pi / 180
theta3 <- 25 * pi / 180

beta_com2 <- cos(theta2) * v1_unit +
  sin(theta2) * v2_unit

beta_com3 <- cos(theta3) * v1_unit -
  sin(theta3) * v2_unit

# Orthogonal projections onto v1
projection2_length <- sum(beta_com2 * v1_unit)
projection3_length <- sum(beta_com3 * v1_unit)

projection2 <- projection2_length * v1_unit
projection3 <- projection3_length * v1_unit

# Colours
com2_colour <- "#009E73"  # green
com3_colour <- "#8E44AD"  # purple

# Angle arcs
make_arc <- function(angle_values, radius) {
  data.frame(
    x = radius * (
      cos(angle_values) * v1_unit[1] +
        sin(angle_values) * v2_unit[1]
    ),
    y = radius * (
      cos(angle_values) * v1_unit[2] +
        sin(angle_values) * v2_unit[2]
    )
  )
}

arc2_df <- make_arc(
  seq(0, theta2, length.out = 100),
  radius = 0.34
)

arc3_df <- make_arc(
  seq(0, -theta3, length.out = 100),
  radius = 0.48
)

theta2_position <- make_arc(theta2 / 2, radius = 0.43)
theta3_position <- make_arc(-theta3 / 2, radius = 0.57)

# Positions of projection-length labels
projection2_label <- 0.60 * projection2 -
  0.25 * v2_unit

projection3_label <- 0.78 * projection3 +
  0.18 * v2_unit

# Plot
main_plot <- ggplot(
  df,
  aes(x = beta1, y = beta2, z = value)
) +
  
  # Quadratic-form contours
  geom_contour(
    colour = "black",
    bins = 20,
    linewidth = 0.45
  ) +
  
  # Dashed unit circle
  geom_path(
    data = circle_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = "blue",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  
  # Projection lengths along v1
  # Purple is drawn first so the shorter green projection remains visible
  annotate(
    "segment",
    x = 0, y = 0,
    xend = projection3[1],
    yend = projection3[2],
    colour = com3_colour,
    linewidth = 5,
    alpha = 0.45
  ) +
  annotate(
    "segment",
    x = 0, y = 0,
    xend = projection2[1],
    yend = projection2[2],
    colour = com2_colour,
    linewidth = 5,
    alpha = 0.75
  ) +
  
  # Perpendicular projection lines
  annotate(
    "segment",
    x = beta_com2[1],
    y = beta_com2[2],
    xend = projection2[1],
    yend = projection2[2],
    colour = com2_colour,
    linetype = "dotted",
    linewidth = 1
  ) +
  annotate(
    "segment",
    x = beta_com3[1],
    y = beta_com3[2],
    xend = projection3[1],
    yend = projection3[2],
    colour = com3_colour,
    linetype = "dotted",
    linewidth = 1
  ) +
  
  # Angle arcs
  geom_path(
    data = arc2_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = com2_colour,
    linewidth = 1
  ) +
  geom_path(
    data = arc3_df,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    colour = com3_colour,
    linewidth = 1
  ) +
  
  # Eigenvectors v1 and v2
  annotate(
    "segment",
    x = 0, y = 0,
    xend = v1[1], yend = v1[2],
    arrow = arrow(
      length = unit(0.28, "cm"),
      type = "closed"
    ),
    colour = "red",
    linewidth = 1
  ) +
  annotate(
    "segment",
    x = 0, y = 0,
    xend = v2[1], yend = v2[2],
    arrow = arrow(
      length = unit(0.28, "cm"),
      type = "closed"
    ),
    colour = "red",
    linewidth = 1
  ) +
  
  # Community vectors
  annotate(
    "segment",
    x = 0, y = 0,
    xend = beta_com2[1],
    yend = beta_com2[2],
    arrow = arrow(
      length = unit(0.25, "cm"),
      type = "closed"
    ),
    colour = com2_colour,
    linewidth = 1.25
  ) +
  annotate(
    "segment",
    x = 0, y = 0,
    xend = beta_com3[1],
    yend = beta_com3[2],
    arrow = arrow(
      length = unit(0.25, "cm"),
      type = "closed"
    ),
    colour = com3_colour,
    linewidth = 1.25
  ) +
  
  # Points at the origin and community-vector endpoints
  annotate(
    "point",
    x = 0, y = 0,
    size = 3.2,
    colour = "black"
  ) +
  annotate(
    "point",
    x = beta_com2[1],
    y = beta_com2[2],
    size = 3,
    colour = com2_colour
  ) +
  annotate(
    "point",
    x = beta_com3[1],
    y = beta_com3[2],
    size = 3,
    colour = com3_colour
  ) +
  
  # Community labels
  annotate(
    "text",
    x = -0.16, y = -0.12,
    label = "com 1",
    colour = "black",
    size = 4.5
  ) +
  annotate(
    "text",
    x = 1.13 * beta_com2[1],
    y = 1.13 * beta_com2[2],
    label = "com 2",
    colour = com2_colour,
    fontface = "bold",
    size = 4.5
  ) +
  annotate(
    "text",
    x = 1.14 * beta_com3[1],
    y = 1.14 * beta_com3[2],
    label = "com 3",
    colour = com3_colour,
    fontface = "bold",
    size = 4.5
  ) +
  
  # Eigenvector labels
  annotate(
    "text",
    x = 1.12 * v1[1],
    y = 1.12 * v1[2],
    label = "bold(v)[1]",
    parse = TRUE,
    colour = "red",
    size = 5
  ) +
  annotate(
    "text",
    x = 1.12 * v2[1],
    y = 1.12 * v2[2],
    label = "bold(v)[2]",
    parse = TRUE,
    colour = "red",
    size = 5
  ) +
  
  # Angle labels
  annotate(
    "text",
    x = theta2_position$x,
    y = theta2_position$y,
    label = "theta[2]",
    parse = TRUE,
    colour = com2_colour,
    size = 4.5
  ) +
  annotate(
    "text",
    x = theta3_position$x,
    y = theta3_position$y,
    label = "theta[3]",
    parse = TRUE,
    colour = com3_colour,
    size = 4.5
  ) +
  
  # Projection-length labels
  annotate(
    "text",
    x = projection2_label[1],
    y = projection2_label[2],
    label =
      'group("|", beta["com 2"] * cos(theta[2]), "|")',
    parse = TRUE,
    colour = com2_colour,
    size = 4
  ) +
  annotate(
    "text",
    x = projection3_label[1],
    y = projection3_label[2],
    label =
      'group("|", beta["com 3"] * cos(theta[3]), "|")',
    parse = TRUE,
    colour = com3_colour,
    size = 4
  ) +
  
  labs(
    x = expression(beta[1]),
    y = expression(beta[2]),
    title = "Quadratic-form contours and community projections"
  ) +
  
  coord_fixed(
    xlim = c(-2, 2),
    ylim = c(-2, 2),
    clip = "off"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank()
  )

print(main_plot)
