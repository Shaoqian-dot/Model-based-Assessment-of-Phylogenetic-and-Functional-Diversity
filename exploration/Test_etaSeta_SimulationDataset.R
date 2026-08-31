# =============================================================================
# Model-based assessment of phylogenetic diversity
#
# This script:
#   1. Generates the simulated community data.
#   2. Fits a reduced-rank GLMM.
#   3. Reconstructs model-based species-composition effects.
#   4. Projects these effects onto phylogenetic eigenvectors.
#   5. Compares conventional Rao's Q with model-based diversity measures.
# =============================================================================
rm(list = ls())

# =============================================================================
# 0. Packages
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
    abundance, # abundance is a matrix which has no community labels.
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

# -----------------------------------------------------------------------------
# 1.1 Calculate eta^T S eta for each observation
# -----------------------------------------------------------------------------
#
# Arguments:
#   abundance : matrix of observations x species
#   S         : species-by-species similarity matrix
#
# Returns:
#   A vector containing eta^T S eta for each observation.
#
# The species order in abundance must match the row/column order of S.

calculate_eta_S_eta <- function(abundance, S) {
  
  # Check dimensions.
  if (
    ncol(abundance) != nrow(S) ||
    ncol(abundance) != ncol(S)
  ) {
    stop(
      "The number of species in abundance must match the dimensions of S."
    )
  }
  
  # Check species names if available.
  if (
    !is.null(colnames(abundance)) &&
    !is.null(rownames(S))
  ) {
    if (!all(colnames(abundance) == rownames(S))) {
      stop(
        "Species order in abundance and S does not match."
      )
    }
  }
  
  # Calculate eta^T S eta for every row.
  abundance <- as.matrix(abundance)
  storage.mode(abundance) <- "numeric"
  result <- rowSums(
    (abundance %*% S) * abundance
  )
  
  return(result)
}


# -----------------------------------------------------------------------------
# 1.2 Boxplot with Dunn's test and compact letter display
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


# -----------------------------------------------------------------------------
# 1.3 Biplot of model-based species-composition scores
# -----------------------------------------------------------------------------
#
# Arguments:
#   abundance : matrix whose first column is "com" and remaining columns
#               contain species-level effects
#   P         : eigenvector matrix
#
# The species-level effects are projected onto the first two eigenvectors.

Biplot <- function(abundance, P) {
  
  # ---------------------------------------------------------------------------
  # Extract species-level effects
  # ---------------------------------------------------------------------------
  
  eta_M <- as.matrix(
    abundance[, setdiff(colnames(abundance), "com")]
  )
  
  # Ensure that the matrix is numeric.
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
    geom_point(size = 3) +
    labs(
      x = expression(v[1]^T * eta),
      y = expression(v[2]^T * eta)
    ) +
    coord_fixed(ratio = 1) +
    theme_classic()
}
# -----------------------------------------------------------------------------
# 1.4 Simulate_random_effects
# -----------------------------------------------------------------------------
simulate_random_effects <- function(model) {
  
  # Conditional means of random effects
  ranef_model <- ranef(model, condVar = TRUE)
  re_mean <- ranef_model$cond$id
  
  # Conditional variance-covariance matrices
  re_var <- attr(re_mean, "condVar")
  
  # Number of IDs and random-effect dimensions
  n <- nrow(re_mean)
  p <- ncol(re_mean)
  
  # Matrix to store simulated random effects
  RR <- matrix(
    NA_real_,
    nrow = n,
    ncol = p
  )
  
  # Generate one multivariate-normal draw for each ID
  for (j in seq_len(n)) {
    
    # Conditional mean for ID j
    mu_j <- as.numeric(re_mean[j, ])
    
    # Conditional covariance for ID j
    Sigma_j <- re_var[, , j]
    Sigma_j[is.na(Sigma_j)] <- 0
    
    # Generate one random-effect vector
    RR[j, ] <- MASS::mvrnorm(
      n = 1,
      mu = mu_j,
      Sigma = Sigma_j
    )
  }
  
  # Keep random-effect names
  colnames(RR) <- colnames(re_mean)
  
  return(RR)
}


# Create a function for producing a mean matrix
# U() must be defined externally and return a matrix of dimension (q*r) x p
mean_matrix <- function (Index_change, beta, alpha, 
                         q, r, NOPS, Corr, P, Eigen, Distribution, p){
  
  # Construct (q-1) rows, each with p elements, alternating between -beta and beta
  Means_logit_beta <- rep(rep(c(-beta, beta), length.out = p), times = q - 1)
  
  # Construct the last row (community 5), each with p elements,
  # alternating between (alpha - beta) and (alpha + beta)
  Means_logit_alpha <- rep(c(alpha - beta, alpha + beta), length.out = p)
  
  # Concatenate all rows into a single vector
  Means_logit <- c(Means_logit_beta, Means_logit_alpha)
  
  # Reshape into a q x p matrix (row-wise filling)
  Mean_matrix_logit <- matrix(Means_logit, c(q, p), byrow = TRUE)
  
  # Loop over groups (excluding the first and last group)
  for (i in 2 : (q - 1)){
    
    # Number of available index pairs for swapping in group (i-1)
    available_pairs <- length(Index_change[[i-1]]) / 2
    
    # If no valid pairs exist, print a message
    if(available_pairs == 0){
      message(sprintf("No satisfied pair in group %d", i-1))
    } else {
      # Number of pairs to sample (cannot exceed available pairs)
      n_sample <- min(NOPS, available_pairs)
      
      # Sample index pairs from the provided index set
      Pairs <- sample_pairs(vec = Index_change[[i-1]], n_sample = n_sample)
      
      # Flip the sign of selected positions in the mean matrix
      Mean_matrix_logit[i, Pairs] <- -Mean_matrix_logit[i, Pairs]
    }
  }
  
  # If correlation structure is required
  if (Corr == 1){
    
    # Generate random effect matrix U_M of dimension (q*r) x p
    U_M <- U(p = p, r = r, q = q, P = P, Eigen = Eigen)
    
    # Expand pean patrix to (q*r) x p and add random effects
    # Then apply inverse logit transformation
    Mean_matrix_logit_cor <- Mean_matrix_logit[rep(1 : q, each = r), ] + U_M
    
  } else {
    # No correlation: only expand and apply inverse link function
    Mean_matrix_logit_cor <- Mean_matrix_logit[rep(1 : q, each = r), ]
  }
  
  if (Distribution == 'binomial') {
    Mean_matrix_full <- plogis(Mean_matrix_logit_cor)
  } else if (Distribution == 'poisson'){
    Mean_matrix_full <- exp(Mean_matrix_logit_cor)
  }
  assign('Mean_matrix_full', Mean_matrix_full, envir =  .GlobalEnv)
  assign('Mean_matrix_logit', Mean_matrix_logit, envir =  .GlobalEnv)
  assign('Mean_matrix_logit_cor', Mean_matrix_logit_cor, envir =  .GlobalEnv)
  # Return final mean matrix
  return(Mean_matrix_full)
}

# -------------------------------------------------------------------------
# 11. The following function chooses the first n_sample pairs of vec instead of 
# randomly sampling the n_sample pairs
# -------------------------------------------------------------------------
# vec: the original vector (length must be even)
# n_sample: number of pairs to sample

sample_pairs <- function(vec, n_sample){
  # Indices of the first element in each pair (1, 3, 5, ...)
  odd_numbers <- seq(1, length(vec), by = 2)                
  
  # Total number of pairs
  n_pairs <- length(vec) / 2 
  
  # Randomly select indices corresponding to pairs
  #selected_pair_idx <- sample(odd_numbers, n_sample)
  selected_pair_idx <- odd_numbers[seq_len(n_sample)]
  # Get the indices of both elements in each selected pair
  selected_elements_idx <- as.vector(
    sapply(selected_pair_idx, function(i) c(i, i + 1))
  )
  
  # Return the selected elements
  return(vec[selected_elements_idx])
}

# -----------------------------------------------------------------------------
# 1.5 test \mu^{\hat}S\mu^{\hat}, \mu^{\hat} = xB^{\hat}
# -----------------------------------------------------------------------------
test_mu_pvalues <- function(p, q, r, NOPS, tree, alpha, beta, quantile, 
                            Corr, Eigen, family_type, 
                            nsim, mu_type = c("mu", "mu_B")) {
  mu_type <- match.arg(mu_type)
  # ---------------------------------------------------------------------------
  # 6. Storage for p-values
  # ---------------------------------------------------------------------------
  
  # First set of tests:
  #   mu^T S mu
  #
  # Second set of tests:
  #   sum_{i=1}^k lambda_i (v_i^T mu)^2
  
  p_values_S <-  rep(NA, q - 1)
  
  col_names <- paste0(
    "com",
    seq_len(q)[-1],
    "_vs_com1"
  )
  
  names(p_values_S) <- col_names
  
  # ---------------------------------------------------------------------------
  # 7. Monte Carlo simulation
  # ---------------------------------------------------------------------------
  DM_phy_func <- getPhyloMatrix(
    tree = tree,
    m = p
  )
  # =============================================================================
  # 5. Construct phylogenetic distance and similarity matrices
  # =============================================================================
  
  # ---------------------------------------------------------------------------
  # Phylogenetic distance matrix
  # ---------------------------------------------------------------------------
  
  
  
  colnames(DM_phy_func) <- paste0(
    "sp",
    seq_len(p)
  )
  
  rownames(DM_phy_func) <- paste0(
    "sp",
    seq_len(p)
  )
  
  # ---------------------------------------------------------------------------
  # Convert phylogenetic distance to a similarity-like matrix
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
  # 6. Centre the phylogenetic similarity matrix
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
  
  test_data <- data.frame(com = paste0(1:5)) # Data frame for t-tests
  for (b in seq_len(nsim)) {
    # =============================================================================
    # 7. Generate community data
    # =============================================================================
    
    my.sample.filtered <- getData(
      DM_phy_func = DM_phy_func,
      NOPS = NOPS,
      q = q,
      alpha = alpha,
      beta = beta,
      r = r,
      quantile = quantile,
      Corr = Corr,
      P = P_J[, 1:(p - 1)],
      Eigen = Eigen,
      Distribution = family_type,
      p = p
    )
    
    # =============================================================================
    # 8. Prepare data for model fitting
    # =============================================================================
    
    yX <- my.sample.filtered$yX
    
    # Convert categorical variables to factors.
    yX[c("sp", "com", "id")] <-
      lapply(
        yX[c("sp", "com", "id")],
        factor
      )
    # ---------------------------------------------------------------------------
    # 4. Community information
    # ---------------------------------------------------------------------------
    com_levels <- levels(yX$com)
    # =============================================================================
    # 9. Fit the reduced-rank GLMM
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
      # control = glmmTMBControl(
      #   optCtrl = list(
      #     iter.max = 1000,
      #     eval.max = 1000
      #   )
      # ),
      data = yX
    )
    
    S <- VC_phy_func # Specify VC_phy_func as S in muSmu.
    
    # ---------------------------------------------------------------------------
    # 2. Fixed-effect predictions
    # ---------------------------------------------------------------------------
    
    eta_FE <- predict(
      fit_Model_5,
      type = "link",
      re.form = NA
    )
    
    eta_FE <- matrix(
      eta_FE,
      ncol = p,
      byrow = FALSE
    )
    
    # Check dimensions
    if (ncol(eta_FE) != ncol(S)) {
      stop(
        "The number of columns in fixed-effect predictions must ",
        "match the dimensions of S."
      )
    }
    
    if (nrow(S) != ncol(S)) {
      stop("S must be a square matrix.")
    }
    
    # ---------------------------------------------------------------------------
    # 3. Row-centred fixed effects for mu_B
    # ---------------------------------------------------------------------------
    
    eta_FE_B <- eta_FE - rowMeans(eta_FE)
    
    # -------------------------------------------------------------------------
    # Construct mu or mu_B
    # -------------------------------------------------------------------------
    
    if (mu_type == "mu") {
      mu <- eta_FE 
    } else {
      mu <- eta_FE_B
    }
    mu <- as.matrix(mu)
    mu <- unique(mu)
    # -------------------------------------------------------------------------
    # 8. Calculate eta^T S eta
    # -------------------------------------------------------------------------
    
    mu_S_mu <- rowSums(
      (mu %*% S) * mu
    )
    
    test_data[[paste0('muSmu', b)]] <- mu_S_mu
  }
  # -------------------------------------------------------------------------
  # 11. t-tests: com 1 versus every other community
  # -------------------------------------------------------------------------
  
  for (j in 2:q) {
    comparison_data <- subset(
      test_data,
      com %in% c(com_levels[1], com_levels[j])
    )
    
    # Test using mu^T S mu
    p_values_S[j - 1] <- t.test(
      comparison_data[1, -1], comparison_data[2, -1]
    )$p.value
  }
  
  assign('eta_FE', eta_FE, envir =  .GlobalEnv)
  assign('eta_FE_B', eta_FE_B, envir =  .GlobalEnv)
  assign('mu',mu, envir =  .GlobalEnv)
  assign('test_data', test_data, envir =  .GlobalEnv)
  assign('VC_phy_func', VC_phy_func, envir =  .GlobalEnv)
  return(p_values_S)
}

# =============================================================================
# 3. File paths and source functions
# =============================================================================

source('R/spectral_decomp.R')
source('R/getData_current.R')
source('R/GroupingPairs.R')
#source('R/mean_matrix_current.R')
#source('R/sample_pairs.R')
source('R/U_current.R')
source('R/getPhyloMatrix.R')
source("R/get_num.eig.R") 
# =============================================================================
# 4. Simulation settings
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

p_values <- test_mu_pvalues(p = p, q = q, r = r, NOPS = NOPS, tree = tree, 
                       alpha = alpha, beta = beta, quantile = quantile, 
          Corr = Corr, Eigen = Eigen, family_type = family_type, 
          nsim = 50, mu_type = c("mu_B"))

matplot(
  t(test_data[, -1]),
  type = "l",
  lty = 1,
  lwd = 2,
  xlab = "Column",
  ylab = expression(mu[B]^T * S * mu[B]),
  ylim  = c(150, 1050)
)
legend(
  "topright",
  legend = test_data$com,
  col = seq_along(test_data$com),
  lty = 1
)
## True muSmu
rowSums(Mean_matrix_logit %*% VC_phy_func * Mean_matrix_logit)

# ######################################################## t.test for real sample
# com <- rep(levels(yX$com), each = r)
# eta_real <- cbind(com, Mean_matrix_logit_cor)
# storage.mode(eta_real) <- 'numeric'
# eta_real[com == 5, 2 : (p + 1)] <- eta_real[com == 5, 2 : (p + 1)] - alpha
# etaSeta_real <- rowSums(eta_real[, 2 : (p + 1)] %*% VC_phy_func_J *  eta_real[, 2 : (p + 1)])
# etaSeta_real_data <- data.frame(
#   com = com,
#   etaSeta = etaSeta_real
# )
# boxplot(etaSeta~com, data = etaSeta_real_data)
# t.test(
#   etaSeta ~ com,
#   data = subset(etaSeta_real_data, com %in% c("1", "4"))
# )
# 
