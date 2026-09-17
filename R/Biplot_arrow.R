Biplot <- function(
    abundance,
    P,
    plot_labs = labs(
      x = expression(v[1]^T * eta[B]),
      y = expression(v[2]^T * eta[B])
    ),
    c = 3 # A multiplier for scaling factor loadings so that the scale of them is compatible with priciple components 
) {
  
  # ---------------------------------------------------------------------------
  # Extract species-level effects
  # ---------------------------------------------------------------------------
  
  eta_M <- as.matrix(
    abundance[, setdiff(colnames(abundance), "com")]
  )
  
  # Ensure that the matrix is numeric
  eta_M <- matrix(
    as.numeric(eta_M),
    nrow = nrow(eta_M),
    ncol = ncol(eta_M),
    dimnames = dimnames(eta_M)
  )
  
  
  # ---------------------------------------------------------------------------
  # Project species-level effects onto the first two eigenvectors
  # ---------------------------------------------------------------------------
  
  aim <- eta_M %*% P[, c(1, 2)]
  
  colnames(aim) <- paste0(
    "v",
    seq_len(ncol(aim))
  )
  
  
  # ---------------------------------------------------------------------------
  # Create plotting data
  # ---------------------------------------------------------------------------
  
  df_scores <- data.frame(
    com = abundance[, "com"],
    score1 = aim[, 1],
    score2 = aim[, 2]
  )
  
  df_scores$com <- factor(
    df_scores$com
  )
  
  
  # ---------------------------------------------------------------------------
  # Create arrow data from the first two columns of P
  # ---------------------------------------------------------------------------
  
  df_arrows <- data.frame(
    x = 0,
    y = 0,
    xend = c * P[, 1],
    yend = c * P[, 2],
    label = paste('sp', 1:nrow(P))
  )
  
  
  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  
  ggplot(
    df_scores,
    aes(
      x = score1,
      y = score2,
      color = com
    )
  ) +
    
    # Species scores
    geom_point(
      size = 3
    ) +
    
    # Arrows
    geom_segment(
      data = df_arrows,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE,
      arrow = arrow(
        length = unit(0.15, "cm"),
        type = "closed"
      ),
      linewidth = 0.5,
      color = "black"
    ) +
    
    # Species labels at the end of arrows
    geom_text(
      data = df_arrows,
      aes(
        x = xend,
        y = yend,
        label = label
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 4,
      hjust = 0.5,
      vjust = -0.5
    ) +
    
    # Center at (0, 0)
    geom_point(
      data = data.frame(
        x = 0,
        y = 0
      ),
      aes(
        x = x,
        y = y
      ),
      inherit.aes = FALSE,
      shape = 4,
      size = 4,
      stroke = 1.2,
      color = "red"
    ) +
    
    # Center label
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "(0, 0)",
      color = "red",
      hjust = -0.2,
      vjust = -0.8,
      size = 4
    ) +
    
    plot_labs +
    coord_fixed(ratio = 1) +
    theme_classic()
}