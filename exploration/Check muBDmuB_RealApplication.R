## Check what is \hat{\mu}_B^{T} D \hat{\mu}_B?
#eta_B
eta_B_hat_data <- data.frame(
  Community = my.sample.filtered$Community,
  abundnace = eta_B_hat
)  

colnames(eta_B_hat_data)[-1] <- selected_species
num_cols <- sapply(eta_B_hat_data, is.numeric)

#mu_B
mu_eta_B_hat_data <- eta_B_hat_data
mu_eta_B_hat_data[num_cols] <- 
  expit(eta_B_hat_data[num_cols])
mu_eta_B_hat_unique <- unique(mu_eta_B_hat_data)
rownames(mu_eta_B_hat_unique) <- mu_eta_B_hat_unique$Community
# Reorder species to conform the order with the order in distance matrix
mu_eta_B_hat_unique <- as.matrix(mu_eta_B_hat_unique[, colnames(DM_phy_func)])
# Diversity
Div_model <- rowSums(mu_eta_B_hat_unique %*% DM_phy_func * mu_eta_B_hat_unique)
print(Div_model)

## Check \hat{\mu}_B^{T} D \hat{\mu}_B manually.
df_mean_Abun <- my.sample.filtered %>%
  group_by(Community) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    )
  )
# Identify numeric columns
num_cols <- sapply(df_mean_Abun, is.numeric)
df_mean_RAbun <- df_mean_Abun
df_mean_RAbun[num_cols] <- 
  expit(
    as.matrix(logit(
      df_mean_Abun[num_cols] / 100 + 0.00001
    )) %*% J
  )
df_mean_RAbun$Community <- factor(
  df_mean_RAbun$Community,
  levels = c('LF', 'MF', 'SF', 'NE', 'SE', 'CG', 'OP', 'OA')
)
df_mean_RAbun <- df_mean_RAbun[order(df_mean_RAbun$Community), ]
df_mean_RAbun <- df_mean_RAbun[, c("Community", colnames(DM_phy_func))]
mean_RAbun <- as.matrix(df_mean_RAbun[, -1])
rownames(mean_RAbun) <- c('LF', 'MF', 'SF', 'NE', 'SE', 'CG', 'OP', 'OA')
Div_manual <- rowSums(mean_RAbun %*% DM_phy_func * mean_RAbun)
print(Div_manual)
