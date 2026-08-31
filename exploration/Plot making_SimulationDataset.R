# =============================================================================
# Section header
# =============================================================================

# ---------------------------------------------------------------------------
# Subsection header
# ---------------------------------------------------------------------------

rm(list = ls())
# =============================================================================
# 1. Packages and source functions
# =============================================================================

library(glmmTMB)
library(ape)
library(picante)
library(MASS)
library(tidyverse)
library(dplyr)
library(FD)
library(DHARMa)
library(ecostats)
library(vegan)
library(emmeans)
library(rstatix)
library(multcompView)
library(ggplot2)

source('R/spectral_decomp.R')
source('R/getData_current.R')
source('R/GroupingPairs.R')
source('R/mean_matrix_current.R')
source('R/sample_pairs.R')
source('R/U_current.R')
source('R/getPhyloMatrix.R')
source("R/get_num.eig.R") 
source('R/RaoQ.R')
source('R/randomization.R')
source('R/rand.RaoQ.fun.R')
source('R/Boxplot.R')
source('R/Biplot.R')
source('R/calculate_eta_S_eta.R')
source('R/simulate_random_effects.R')
# =============================================================================
# 1. Simulation settings
# =============================================================================
set.seed(1)
tree <- read.tree("data/example.tre")

family_type <- "poisson"

p <- 5
q <- 5
r <- 128

alpha <- log(9 / 4)
beta <- log(2)

quantile <- 0.35
Corr <- 1
Eigen <- 2
NOPS <- 1
DM_ge <- 1 # 1: DM is generated from the PCV data set; 2 DM is manually designed. 

# =============================================================================
# 2. Construct phylogenetic distance and similarity matrices
# =============================================================================

# ---------------------------------------------------------------------------
# Phylogenetic distance matrix and similarity matrix
# ---------------------------------------------------------------------------
species <- paste0(
  "sp",
  seq_len(p)
)
if (DM_ge == 1) {
  DM_phy_func <- getPhyloMatrix(
    tree = tree,
    m = p
  )
} else if (DM_ge == 2){
  DM_phy_func <- matrix(c(0, 1, 2, 3, 4,
                          1, 0, 2, 3, 4,
                          2, 2, 0, 3, 4,
                          3, 3, 3, 0, 4,
                          4, 4, 4, 4, 0), p, p)
}

colnames(DM_phy_func) <- species
rownames(DM_phy_func) <- species


# ---------------------------------------------------------------------------
# Convert phylogenetic distance to a similarity matrix
# ---------------------------------------------------------------------------

VC_phy_func <- max(DM_phy_func) - DM_phy_func


# ---------------------------------------------------------------------------
# Spectral decomposition
# ---------------------------------------------------------------------------

spec_VC_phy_func <- spectral_decomp(
  VC_phy_func = VC_phy_func
)

P <- spec_VC_phy_func$P
D <- spec_VC_phy_func$D


# =============================================================================
# 3. Centre the phylogenetic similarity matrix
# =============================================================================

# The centering matrix is:
#
#   J = I - 11^T / p
#
# This projects the similarity matrix onto the subspace orthogonal to the
# constant species vector. The resulting eigenvectors therefore represent
# contrasts among species.

J <- diag(p) -
  (1 / p) * matrix(1, p, p)

VC_phy_func_J <- J %*%
  VC_phy_func %*%
  J

# ---------------------------------------------------------------------------
# Spectral decomposition of the centred similarity matrix
# ---------------------------------------------------------------------------

spec_VC_phy_func_J <- spectral_decomp(
  VC_phy_func = VC_phy_func_J
)

P_J <- spec_VC_phy_func_J$P
D_J <- spec_VC_phy_func_J$D


# =============================================================================
# 4. Generate community data
# =============================================================================

my.sample <- getData(
  DM_phy_func = DM_phy_func,
  NOPS = NOPS,
  q = q,
  alpha = alpha,
  beta = beta,
  r = r,
  quantile = quantile,
  Corr = Corr,
  P = P_J[, -1],
  Eigen = Eigen,
  Distribution = family_type,
  p = p
)

# =============================================================================
# 5. Prepare data for model fitting
# =============================================================================

yX <- my.sample$yX

# Convert categorical variables to factors.
yX[c("sp", "com", "id")] <-
  lapply(
    yX[c("sp", "com", "id")],
    factor
  )

comLevel <- levels(yX$com) 
comOrder <- c('1', '2', '3', '4', '5') # Put the desired order of all the levels
com <- rep(comLevel, each = r)
# =============================================================================
# 6. Fit the reduced-rank GLMM
# =============================================================================

# Model:
#
#   y ~ species
#       + community
#       + species × community
#       + reduced-rank random effects
#
# The reduced-rank component captures residual covariance among species.

fit_Model_5 <- glmmTMB(
  y ~ sp +
    com +
    sp:com +
    rr(sp + 0 | id, 2),
  family = "poisson",
  REML = FALSE,
  data = yX
)

# # =============================================================================
# # 7. Extract fixed effects
# # =============================================================================


# ---------------------------------------------------------------------------
# Fixed effects with total and relative abundance
# ---------------------------------------------------------------------------

eta_TA_fe <- predict(
  fit_Model_5,
  type = "link",
  re.form = NA
)
eta_TA_fe_M <- matrix(
  as.numeric(eta_TA_fe),
  ncol = p,
  byrow = FALSE
)

# ---------------------------------------------------------------------------
# Fixed effects with relative abundance
# ---------------------------------------------------------------------------
eta_RA_fe_M <- eta_TA_fe_M - matrix(apply(eta_TA_fe_M, 1, mean), ncol = p, nrow = nrow(eta_TA_fe_M), byrow = FALSE)

#colnames(eta_RA_fe_M) <- species

# =============================================================================
# 8. Reduced-rank residual effects
# =============================================================================
RR <- simulate_random_effects(model = fit_Model_5)

# =============================================================================
# 9. Construct model-based species composition
# =============================================================================

# ---------------------------------------------------------------------------
# Obtain \eta_B
# ---------------------------------------------------------------------------

eta_RA_M <- eta_RA_fe_M + RR

# Add community labels as the first column.
eta_RA_M <- cbind(
  com = com,
  eta_RA_M
)

# ---------------------------------------------------------------------------
# Obtain \eta_T
# ---------------------------------------------------------------------------

# Obtain the fitted linear predictor for each observation.
eta_TA_M <- eta_TA_fe_M + RR

# Add community labels.
eta_TA_M <- cbind(
  com = com,
  eta_TA_M
)

# =============================================================================
# 10. Biplot of model-based species composition
# =============================================================================

Biplot(
  abundance = eta_RA_M,
  P = P
)

Biplot(
  abundance = eta_TA_M,
  P = P
)

# =============================================================================
# 11. Model-based diversity: eta_B^T S eta_B and eta_T^T S eta_T
# =============================================================================

# ---------------------------------------------------------------------------
# Calculate eta_B^T S eta_B
# ---------------------------------------------------------------------------

DivModel_etaB <- calculate_eta_S_eta(
  abundance = eta_RA_M,
  S = VC_phy_func
)


# ---------------------------------------------------------------------------
# Calculate eta_T^T S eta_T
# ---------------------------------------------------------------------------
DivModel_etaT <- calculate_eta_S_eta(
  abundance = eta_TA_M,
  S = VC_phy_func
)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig1 <- Boxplot(
  data = DivModel_etaB,
  group_levels = comOrder,
  ylab = expression(eta[B]^T * S * eta[B])
)

fig1
fig2 <- Boxplot(
  data = DivModel_etaT,
  group_levels = comOrder,
  ylab = expression(eta[T]^T * S * eta[T])
)
fig2


# =============================================================================
# 12. Conventional Rao's Q
# =============================================================================

abundance <- my.sample$abundance

# ---------------------------------------------------------------------------
# Randomized Rao's Q
# ---------------------------------------------------------------------------

RaosQ_rand <- randomization(
  abundance = abundance,
  DM_phy_func = DM_phy_func
)

RaosQ_rand_data <- data.frame(
  group = com,
  value = RaosQ_rand
)



# ---------------------------------------------------------------------------
# Rao's Q
# ---------------------------------------------------------------------------

RaosQ <- RaoQ(
  my.sample = abundance,
  DM = DM_phy_func
)

RaosQ_data <- data.frame(
  group = com,
  value = RaosQ
)

# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

fig3 <- Boxplot(
  data = RaosQ_rand_data,
  group_levels = comOrder,
  ylab = "Randomized Rao's Q"
)

fig3
fig4 <- Boxplot(
  data = RaosQ_data,
  group_levels = comOrder,
  ylab = "Rao's Q"
)

fig4

# =============================================================================
# 13. True Diversity
# =============================================================================

# \mu S_J \mu
rowSums(Mean_matrix_logit %*% VC_phy_func_J * Mean_matrix_logit)
# \mu D_J \mu
DM_phy_func_J <- J %*% DM_phy_func %*% J
rowSums(Mean_matrix_logit %*% DM_phy_func_J * Mean_matrix_logit)
# \mu S \mu
rowSums(Mean_matrix_logit %*% VC_phy_func * Mean_matrix_logit)
# \mu D \mu
rowSums(Mean_matrix_logit %*% DM_phy_func * Mean_matrix_logit)
# RaosQ
RaoQ(my.sample = exp(Mean_matrix_logit), DM =DM_phy_func)
# Rand RaosQ
randomization(abundance = exp(Mean_matrix_logit), DM_phy_func = DM_phy_func)

