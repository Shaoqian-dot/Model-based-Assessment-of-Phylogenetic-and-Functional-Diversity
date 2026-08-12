###############################################################################
# Model-based Assessment of Phylogenetic and Functional Diversity
#
# Real-data analysis:
#   1. Prepare community and phylogenetic data
#   2. Fit a reduced-rank GLMM
#   3. Extract species-level relative-abundance effects
#   4. Reconstruct model-based species composition
#   5. Project composition onto eigenvectors of the similarity matrix
#   6. Assess group differences using randomized Rao's Q
###############################################################################


# =============================================================================
# 0. Packages
# =============================================================================

library(glmmTMB)
library(ape)
library(picante)
library(MASS)
library(tidyverse)
library(FD)
library(DHARMa)
library(ecostats)
library(vegan)
library(emmeans)
library(rstatix)
library(multcompView)


# =============================================================================
# 1. Helper functions
# =============================================================================

# -----------------------------------------------------------------------------
# Rao's Quadratic Entropy
#
# For each community/site:
#
#   Q_i = 0.5 * p_i' D p_i
#
# where p_i is the vector of relative species abundances and D is a
# species-by-species distance matrix.
# -----------------------------------------------------------------------------

RaoQ <- function(my.sample, DM) {
  
  # Convert abundances to relative abundances within each site.
  relative_abundance <- my.sample / rowSums(my.sample)
  
  # Sites with zero total abundance have undefined relative abundances.
  # Treat these as having zero diversity.
  relative_abundance[is.na(relative_abundance)] <- 0
  
  P <- as.matrix(relative_abundance)
  
  # Vectorized calculation of Rao's Q for all sites.
  RaoQ_diversity <- 0.5 * rowSums((P %*% DM) * P)
  
  return(RaoQ_diversity)
}


# -----------------------------------------------------------------------------
# Randomized Rao's Q
#
# Randomizes species labels on the distance matrix while keeping the
# community abundance matrix unchanged.
#
# This produces one value of Rao's Q under a randomized species-distance
# association.
# -----------------------------------------------------------------------------

rand.RaoQ.fun <- function(DM, my.sample) {
  
  # Species represented in the distance matrix.
  species <- colnames(DM)
  
  # Randomly permute species labels.
  species_permutation <- sample(
    species,
    size = length(species),
    replace = FALSE
  )
  
  # Apply the same permutation to rows and columns.
  DM_shuffle <- DM
  rownames(DM_shuffle) <- species_permutation
  colnames(DM_shuffle) <- species_permutation
  
  # Reorder the abundance matrix to match the randomized distance matrix.
  abundance_shuffle <- my.sample[, colnames(DM_shuffle), drop = FALSE]
  
  # Calculate Rao's Q under the randomized association.
  RaoQ(
    my.sample = abundance_shuffle,
    DM = DM_shuffle
  )
}


# -----------------------------------------------------------------------------
# Randomization test for Rao's Q
#
# Returns the standardized effect size (SES):
#
#   SES = (observed Rao's Q - mean(null Rao's Q)) /
#         sd(null Rao's Q)
#
# The default is 999 randomizations.
# -----------------------------------------------------------------------------

randomization <- function(
    abundance,
    DM_phy_func,
    nperm = 999
) {
  
  # Ensure that abundance columns correspond to the distance matrix.
  colnames(abundance) <- colnames(DM_phy_func)
  
  # Generate the null distribution.
  null_output <- replicate(
    nperm,
    rand.RaoQ.fun(
      DM = DM_phy_func,
      my.sample = abundance
    ),
    simplify = "matrix"
  )
  
  # Observed Rao's Q.
  observed_RaoQ <- RaoQ(
    my.sample = abundance,
    DM = DM_phy_func
  )
  
  # Standardized effect size.
  null_mean <- apply(
    null_output,
    MARGIN = 1,
    mean,
    na.rm = TRUE
  )
  
  null_sd <- apply(
    null_output,
    MARGIN = 1,
    sd,
    na.rm = TRUE
  )
  
  ses <- (observed_RaoQ - null_mean) / null_sd
  
  # If the null distribution has zero variance, SES is undefined.
  # Here it is set to zero.
  ses[is.nan(ses)] <- 0
  
  return(ses)
}


# =============================================================================
# 2. File paths and source functions
# =============================================================================

# NOTE:
# Consider replacing these absolute paths with project-relative paths if
# this analysis will eventually be shared or moved to another computer.

functions_dir <- paste0(
  "C:/Users/huang/OneDrive - UNSW/PhD Project/",
  "Model contrasts/Simulation/Simulation with wide settings/",
  "Code refactoring/Model-based assessment of PD and FD/",
  "Katana/Assess_Diversity_Corr_FL_P_current_bino_pois/functions"
)

data_dir <- paste0(
  "C:/Users/huang/OneDrive - UNSW/PhD Project/",
  "Model contrasts/Simulation/Simulation with wide settings/",
  "Code refactoring/Model-based assessment of PD and FD/",
  "Real application/data"
)

source(file.path(functions_dir, "get_num.eig.R"))
source(file.path(functions_dir, "spectral_decomp.R"))

setwd(data_dir)


# =============================================================================
# 3. Read raw data
# =============================================================================

tree <- read.tree("example.tre")

trait <- read.csv(
  "trait data.csv",
  stringsAsFactors = FALSE
)

my.sample <- read.csv(
  "community data.csv",
  stringsAsFactors = FALSE
)


# =============================================================================
# 4. Prepare trait data
# =============================================================================

# Log-transform and standardize continuous traits.
# The log transformation reduces skewness before standardization.

trait$height_mean <- as.numeric(
  scale(log(trait$height_mean))
)

trait$X1000.S <- as.numeric(
  scale(log(trait$X1000.S))
)

trait$SLA_mean <- as.numeric(
  scale(log(trait$SLA_mean))
)

trait$flowering.duration <- as.numeric(
  scale(log(trait$flowering.duration))
)


# Remove variables identified as highly collinear.
cols_to_drop <- c(
  "flowering_start.7.9",
  "LF_E",
  "seed.disp._ZO",
  "pollination_SELF",
  "reprod_veg"
)

trait <- trait %>%
  select(-all_of(cols_to_drop))


# Replace missing values in trait columns with zero.
# NOTE: This assumes that zero is an appropriate representation of missing
# values for these variables. This should be justified in the Methods section.
trait[, 7:19][is.na(trait[, 7:19])] <- 0


# Convert categorical trait variables to ordered factors.
for (j in 7:19) {
  trait[, j] <- ordered(trait[, j])
}


# =============================================================================
# 5. Harmonize species names between tree and community data
# =============================================================================

# Convert underscores in tree tip labels to periods to match the community data.
tree_species <- str_replace_all(
  tree$tip.label,
  "_",
  "."
)

community_species <- colnames(my.sample)

# Identify tree species whose names do not have a unique match in the
# community dataset.
unmatched_index <- which(
  !tree_species %in% community_species
)

# Manually correct known species-name mismatches.
tree_species[unmatched_index] <- c(
  "Astragalus.dasyanthus.exscapus",
  "Populus.alba...P..x.canescens",
  "Capsella.bursa.pastoris",
  "Festuca.rupicola.valesiaca",
  "Stipa.borysthenica.capillata"
)

tree$tip.label <- tree_species


# =============================================================================
# 6. Filter community data
# =============================================================================

# Remove species occurring in fewer than three plots.
rare_species <- c(
  "Acer.sp.",
  "Epipactis.sp.",
  "Fraxinus.sp.",
  "Hieracium.sp.",
  "Lathyrus.sp.",
  "Prunus.sp.",
  "Silene.sp."
)

my.sample <- my.sample %>%
  select(-any_of(rare_species))


# Calculate basic diversity summaries for exploratory purposes.
SR <- specnumber(
  my.sample[, -c(1, 2)]
)

Shannon_D <- diversity(
  my.sample[, -c(1, 2)],
  index = "shannon"
)


# Remove species occurring in more than 95% of plots.
# This retains species with zero values in no more than 95% of observations.
zero_prop <- colMeans(
  my.sample == 0,
  na.rm = TRUE
)

my.sample.filtered <- my.sample %>%
  select(
    where(~ TRUE)
  )

# Apply the occurrence filter only to species columns.
species_columns <- colnames(my.sample)[-c(1, 2)]

species_to_keep <- species_columns[
  zero_prop[species_columns] <= 0.95
]

my.sample.filtered <- my.sample %>%
  select(
    all_of(colnames(my.sample)[1:2]),
    all_of(species_to_keep)
  )


# Species retained for the model.
selected_species <- colnames(
  my.sample.filtered
)[-c(1, 2)]

m <- length(selected_species)


# =============================================================================
# 7. Convert community data to long format
# =============================================================================

yX <- reshape(
  my.sample.filtered,
  direction = "long",
  varying = selected_species,
  v.names = "y",
  timevar = "sp"
)

# Convert grouping variables to factors.
yX[c("sp", "Site", "Community", "id")] <-
  lapply(
    yX[c("sp", "Site", "Community", "id")],
    factor
  )


# Convert percentage cover to proportions.
yX$y <- yX$y / 100

# ordBeta requires values strictly between 0 and 1.
yX$y[yX$y >= 1] <- 0.9999


# =============================================================================
# 8. Construct phylogenetic distance and similarity matrices
# =============================================================================

# Keep only species present in the filtered community dataset.
subtree <- keep.tip(
  tree,
  selected_species
)


# Phylogenetic distance matrix.
DM_phy_func <- cophenetic.phylo(subtree)


# Convert distance to a similarity-like matrix.
VC_phy_func <- max(DM_phy_func) - DM_phy_func


# Spectral decomposition of the similarity matrix.
spec_VC_phy_func <- spectral_decomp(
  VC_phy_func = VC_phy_func
)

P <- spec_VC_phy_func$P
D <- spec_VC_phy_func$D


# Determine the number of eigenvectors using Method 3.
val_num.eig <- get_num.eig(
  "Method3",
  NA,
  D
)


# =============================================================================
# 9. Project the similarity matrix onto the centred subspace
# =============================================================================

# J is the centering matrix:
#
#   J = I - 11'/m
#
# This removes the constant species component and ensures that the
# eigenvectors describe contrasts among species.

J <- diag(m) -
  (1 / m) * matrix(1, m, m)


VC_phy_func_J <- J %*%
  VC_phy_func %*%
  J


spec_VC_phy_func_J <- spectral_decomp(
  VC_phy_func = VC_phy_func_J
)

P_J <- spec_VC_phy_func_J$P
D_J <- spec_VC_phy_func_J$D


# Determine the number of eigenvectors after centering.
val_num.eig_J <- get_num.eig(
  "Method3",
  NA,
  D_J[1:(m - 1), 1:(m - 1)]
)


# =============================================================================
# 10. Fit the reduced-rank GLMM
# =============================================================================

# Model structure:
#
#   y ~ species + community + species:community
#       + reduced-rank random effects
#
# The reduced-rank component captures residual covariance among species.

fit_Model_5 <- glmmTMB(
  y ~ sp +
    Community +
    sp:Community +
    rr(sp + 0 | id, 2),
  family = "ordbeta",
  REML = FALSE,
  control = glmmTMBControl(
    optCtrl = list(
      iter.max = 1000,
      eval.max = 1000
    )
  ),
  data = yX
)


# =============================================================================
# 11. Extract fixed effects
# =============================================================================

beta <- fixef(
  fit_Model_5
)$cond

species <- levels(yX$sp)
communities <- levels(yX$Community)

m <- length(species)
ncom <- length(communities)


# =============================================================================
# 12. Construct species-level relative-abundance effects
# =============================================================================

# eta_RA contains the model-estimated species effects for each community.
#
# Rows    = communities
# Columns = species
#
# The first species and first community act as reference levels under the
# treatment coding used by the model.

eta_RA <- matrix(
  0,
  nrow = ncom,
  ncol = m,
  dimnames = list(
    communities,
    species
  )
)


# -----------------------------------------------------------------------------
# 12.1 Main species effects
# -----------------------------------------------------------------------------

main_names <- paste0(
  "sp",
  species[-1]
)

# First community is the reference community.
eta_RA[1, 1] <- 0

for (j in 2:m) {
  
  if (main_names[j - 1] %in% names(beta)) {
    eta_RA[1, j] <- beta[
      main_names[j - 1]
    ]
  }
}


# -----------------------------------------------------------------------------
# 12.2 Species-by-community interaction effects
# -----------------------------------------------------------------------------

for (i in 2:ncom) {
  
  current_community <- communities[i]
  
  # First species is the reference species.
  eta_RA[i, 1] <- 0
  
  for (j in 2:m) {
    
    interaction_name <- paste0(
      "sp",
      species[j],
      ":Community",
      current_community
    )
    
    if (interaction_name %in% names(beta)) {
      eta_RA[i, j] <- beta[
        interaction_name
      ]
    }
  }
}


# =============================================================================
# 13. Centre species effects within each community
# =============================================================================

# The diversity contrast is defined in the subspace satisfying:
#
#   sum(beta_j) = 0
#
# Therefore, centre the estimated species effects within each community.

eta_RA <- sweep(
  eta_RA,
  MARGIN = 1,
  STATS = rowMeans(eta_RA),
  FUN = "-"
)


# =============================================================================
# 14. Expand community-level effects to individual observations
# =============================================================================

obs_eta <- eta_RA[
  match(
    my.sample.filtered$Community,
    communities
  ),
  ,
  drop = FALSE
]


# =============================================================================
# 15. Generate reduced-rank residual effects
# =============================================================================

# Extract the fitted reduced-rank factor loadings.
model_report <- fit_Model_5$obj$env$report()

Lambda <- model_report$fact_load[[1]]

rank <- ncol(Lambda)


# Generate standard normal latent variables.
Z <- MASS::mvrnorm(
  n = nrow(obs_eta),
  mu = rep(0, rank),
  Sigma = diag(rank)
)


# Convert latent variables into species-level residual effects.
RR <- Z %*% t(Lambda)


# =============================================================================
# 16. Reconstruct model-based species composition
# =============================================================================

# eta_total represents the model-based species composition:
#
#   eta_total = fixed species/community effects
#               + reduced-rank residual effects
#
# Rows    = observations
# Columns = species

eta_total <- obs_eta + RR


# =============================================================================
# 17. Project model-based composition onto phylogenetic eigenvectors
# =============================================================================

# Project each observation's species-composition vector onto the first two
# eigenvectors of the centred similarity matrix.

eta_M <- eta_total

aim <- eta_M %*%
  P_J[, c(1, 2)]

colnames(aim) <- paste0(
  "p",
  1:ncol(aim)
)


# Create plotting dataset.
df_scores <- data.frame(
  com = my.sample.filtered$Community,
  score1 = aim[, 1],
  score2 = aim[, 2]
)


# Set the desired order of communities.
df_scores$com <- factor(
  df_scores$com,
  levels = c(
    "LF",
    "MF",
    "SF",
    "NE",
    "SE",
    "CG",
    "OP",
    "OA"
  )
)


# =============================================================================
# 18. Plot the model-based species-composition scores
# =============================================================================

ggplot(
  df_scores,
  aes(
    x = score1,
    y = score2,
    color = com
  )
) +
  geom_point(size = 3) +
  labs(
    x = expression(v[1]^T * eta),
    y = expression(v[2]^T * eta)
  ) +
  theme_classic()


# =============================================================================
# 19. Randomized Rao's Q
# =============================================================================

# Ensure that the abundance matrix has the same species and ordering as the
# phylogenetic distance matrix.

abundance <- my.sample.filtered[
  ,
  colnames(DM_phy_func),
  drop = FALSE
]


# Calculate SES of Rao's Q using tip randomization.
RaosQ_rand <- randomization(
  abundance = abundance,
  DM_phy_func = DM_phy_func
)


# Combine Rao's Q values with community groups.
RaosQ_rand_data <- data.frame(
  group = df_scores$com,
  value = RaosQ_rand
)


# =============================================================================
# 20. Dunn's test with BH correction
# =============================================================================

RaosQ_rand_data$group <- factor(
  RaosQ_rand_data$group,
  levels = c(
    "LF",
    "MF",
    "SF",
    "NE",
    "SE",
    "CG",
    "OP",
    "OA"
  )
)


stat.test <- RaosQ_rand_data %>%
  dunn_test(
    value ~ group,
    p.adjust.method = "BH"
  )


# =============================================================================
# 21. Compact letter display
# =============================================================================

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
  levels = levels(RaosQ_rand_data$group)
)


# Position the significance letters above the boxplots.
y_pos <- max(
  RaosQ_rand_data$value,
  na.rm = TRUE
) * 1.1


# =============================================================================
# 22. Plot randomized Rao's Q
# =============================================================================

ggplot(
  RaosQ_rand_data,
  aes(
    x = group,
    y = value,
    fill = group
  )
) +
  
  # Boxplots show the distribution of randomized Rao's Q.
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  
  # Add individual observations.
  geom_jitter(
    width = 0.15,
    alpha = 0.5,
    size = 2
  ) +
  
  # Add compact-letter significance groups.
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
    y = "Randomized Rao's Q"
  )

