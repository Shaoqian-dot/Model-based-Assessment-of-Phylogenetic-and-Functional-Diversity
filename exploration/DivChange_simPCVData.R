rm(list = ls())
library(ape)
library(glmmTMB)
source('R/sample_pairs.R')
source('R/RaoQ.R')
source('R/rand.RaoQ.fun.R')
source('R/randomization.R')
# Function for creating mean matrix with five communities
mean_matrix <- function (beta, alpha, q, NOPS, abund_ref, # abund_ref is the selected species abundances waiting to swap/change
                         swap = FALSE, effectSize = FALSE, tree){ 
  
  p <- length(abund_ref)

  # Alternatively rearrange species abundances into large, small, large, small...  
  abund_ref_order <- sort(abund_ref) # abund_ref is a stratum - subsample from the original 115 species
  
  abund_ref_small <- sample(abund_ref_order[1 : (p/2)])
  abund_ref_large <- sample(abund_ref_order[floor(((p/2) + 1)) : p])
  
  abund_com1 <- interleave(abund_ref_large, abund_ref_small)# The resulting species abudnance presenting as large, small, large, small,...
  species <- names(abund_com1)
  # Swap or change with the same effect size to create abundance matrix
  mean_matrix <- matrix(rep(abund_com1, times = q), 
                        nrow = q, byrow = TRUE)
  mean_matrix_logit <- qlogis(mean_matrix)
  mean_matrix_logit[q, ] <- qlogis(abund_com1) + alpha # Create community 5
  colnames(mean_matrix_logit) <- species
  subtree <- keep.tip(tree, species)
  #plot(subtree)
  ## Compute the pairwise phylogenetic distance matrix (cophenetic distances)
  DM_phy_func <- cophenetic.phylo(subtree)
  ## Generate Index to swap
  Index_change <- GroupingPairs(DM_phy_func = DM_phy_func, quantile = 0.35, species = species)
  
  ## Loop over groups (excluding the first and last group)
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
      
      target <- mean_matrix_logit[i, ]
      if (swap == TRUE) {
        # Flip the mean abundances of selected positions in the mean matrix
        mean_matrix_logit[i, Pairs[seq(1, length(Pairs), by = 2)]] <- target[Pairs[seq(2, length(Pairs), by = 2)]]
        mean_matrix_logit[i, Pairs[seq(2, length(Pairs), by = 2)]] <- target[Pairs[seq(1, length(Pairs), by = 2)]]
      } else if (effectSize == TRUE) {
        mean_matrix_logit[i, Pairs[seq(1, length(Pairs), by = 2)]] <- target[Pairs[seq(1, length(Pairs), by = 2)]] - beta
        mean_matrix_logit[i, Pairs[seq(2, length(Pairs), by = 2)]] <- target[Pairs[seq(1, length(Pairs), by = 2)]] + beta
      }
      
    }
  }
  
  mean_matrix <- plogis(mean_matrix_logit)
  rownames(mean_matrix) <- paste0('com', 1 : q)
  return(mean_matrix)
}

# Function for stratified sampling from a community mean abundance
sample_species_stratified <- function(mean_abundance, p, seed = NULL) {
  set.seed(seed)
  n_species <- length(mean_abundance)
  
  # Order species abundance from rarest to most abundant
  mean_abundance_sort <- sort(mean_abundance)
  
  # Group species abundances 
  strata1 <- mean_abundance_sort[1 : 5]
  strata2 <- mean_abundance_sort[6 : 20]
  strata3 <- mean_abundance_sort[21 : 80]
  strata4 <- mean_abundance_sort[81 : 105]
  strata5 <- mean_abundance_sort[106 : 115]
  
  # Sampling
  abundance_target <- c(sample(strata1, 1), sample(strata2, 1),
                        sample(strata3, 1), sample(strata4, 1),
                        sample(strata5, 1))
  return(abundance_target)
}

# Alternatively arrange x and y
interleave <- function(x, y) {
  
  # Start with the longer vector
  if (length(y) > length(x)) {
    return(interleave(y, x))
  }
  
  result <- c(x, y)
  
  
  result[seq(1, length(result), by = 2)] <- x
  
  result[seq(2, length(result), by = 2)] <- y
  
  species <- names(result)
  species[seq(1, length(result), by = 2)] <- names(x)
  species[seq(2, length(result), by = 2)] <- names(y)
  names(result) <- species
  result
}

# The idea here is to separate 
GroupingPairs <- function(DM_phy_func, quantile, species){# species is the reordered species labels
  # Rerrange the distance matrix according to new species labels
  DM_phy_func_reArrange <- DM_phy_func[species, species]
  
  # Number of species
  m <- length(species)
  
  # Initialize lists (use list for efficiency, convert later if needed)
  Index_No_change   <- integer(0)
  Index_Small_change <- integer(0)
  Index_Big_change   <- integer(0)
  
  # Define groups
  G1 <- seq(1, m, by = 2)
  G2 <- seq(2, m, by = 2)
  
  # Compute distance threshold
  dis_quant <- quantile(DM_phy_func_reArrange[upper.tri(DM_phy_func_reArrange)], quantile)
  
  for (i in G1){
    for (j in G2){
      
      # Case 1: identical phylogenetic profile
      if (isTRUE(all.equal(
        unname(DM_phy_func_reArrange[i, -c(i, j)]),
        unname(DM_phy_func_reArrange[j, -c(i, j)])
      ))) {
        
        if (!any(c(i, j) %in% Index_No_change)) {
          Index_No_change <- c(Index_No_change, i, j)
        }
        
        # Case 2: small change
      } else if (DM_phy_func_reArrange[i, j] <= dis_quant) { # Pairs in small change group exclude those pairs in no change group since if-else if are excluded. 
        
        if (!any(c(i, j) %in% Index_Small_change)) {
          Index_Small_change <- c(Index_Small_change, i, j)
        }
        
        # Case 3: big change
      } else if (DM_phy_func_reArrange[i, j] > dis_quant) { # Pairs in large change group exclude those pairs in no change group since if-else if are excluded. 
        
        if (!any(c(i, j) %in% Index_Big_change)) {
          Index_Big_change <- c(Index_Big_change, i, j)
        }
      }
    }
  }
  
  # Return result (no global assignment!)
  return(list(
    No_change   = Index_No_change,
    Small_change = Index_Small_change,
    Big_change   = Index_Big_change
  ))
}


# 1. Create a new, empty environment
temp_env <- new.env()

# 2. Load your saved RData file into that specific environment
load("results/Real Application/my_environment.RData", envir = temp_env)

# 3. Extract only the specific variable you want into your Global Environment
model <- temp_env$model


model_updated <- glmmTMB::up2date(model)

# 4. Clean up the temporary environment to free memory
rm(temp_env)

alpha <- 2.426285
beta <- 2.426285
q <- 5
NOPS <- 1
Corr <- 0
p_total <- 115
p <- 5 
seed <- 129
swap <- TRUE
effectSize <- FALSE

my.sample <- read.csv(
  "data/community data.csv",
  stringsAsFactors = FALSE)
  

tree <- read.tree("data/example.tre")

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


## eta_T_hat
eta_T_hat <- predict(model_updated, 
                     type = "link",
                     re.form = NA)
eta_T_hat <- matrix(
  as.numeric(eta_T_hat),
  ncol = p_total,
  byrow = FALSE
)
eta_T_hat_unique <- unique(eta_T_hat)

mean_abundance <- plogis(eta_T_hat_unique[8, ])
names(mean_abundance) <- levels(model$frame$sp)

selected_5 <- sample_species_stratified(
  mean_abundance = mean_abundance,
  p = 5,
  seed = seed
)

mu <- mean_matrix(beta = beta, alpha = alpha, q = 5, 
                           NOPS = 1, abund_ref = selected_5, 
                           swap = swap, effectSize = effectSize, tree  = tree)
species <- colnames(mu)
# Distance Matrix
subtree <- keep.tip(tree, species)
plot(subtree)
## Compute the pairwise phylogenetic distance matrix (cophenetic distances)
DM <- cophenetic.phylo(subtree)[species, species]

# Diversity
## Rao's Q
mu_RA <- mu/apply(mu, 1, sum)
Rao <- rowSums(mu_RA %*% DM * mu_RA)
print(Rao)

## Randomized Rao's Q
Rand_RaosQ <- randomization(
  abundance = mu,
  DM_phy_func = DM
)
print(Rand_RaosQ)

## Model mu_B
J <- diag(1, p) - 1/p * matrix(1, p, p)
muB <- plogis(qlogis(mu) %*% J)
Model_muB <- rowSums(muB %*% DM * muB)
print(Model_muB)

if(swap == TRUE & effectSize == FALSE){
  Div_swap <- rbind(Rao, Rand_RaosQ, Model_muB)
} else if(swap == FALSE & effectSize == TRUE) {
  Div_EffectSize <- rbind(Rao, Rand_RaosQ, Model_muB)
}

# Example settings
cols <- c("red", "black")
labs <- c("effectSize", "Swap")

# Three panels in the first row; shared legend in the second row
layout(
  matrix(c(1, 2, 3,
           4, 4, 4),
         nrow = 2, byrow = TRUE),
  heights = c(8, 1)
)

# Plot margins
par(mar = c(4, 4, 2, 1))

# Panel 1
plot(Div_swap[1, ], main = 'Rao', type = 'l', 
     ylim = c(min(c(Div_swap[1, ], Div_EffectSize[1, ])) - 5, 
              max(c(Div_swap[1, ], Div_EffectSize[1, ])) + 5), 
     ylab = 'Diversity', xlab = 'Community')
lines(Div_EffectSize[1, ], col ='red')

# Panel 2
plot(Div_swap[2, ], main = 'Rand', type = 'l',
     ylab = 'Diversity', xlab = 'Community',
     ylim = c(min(c(Div_swap[2, ], Div_EffectSize[2, ])) - 0.5, 
              max(c(Div_swap[2, ], Div_EffectSize[2, ])) + 0.5))
lines(Div_EffectSize[2, ], col = 'red')


# Panel 3
plot(Div_swap[3, ], main = 'Model', type = 'l', 
     ylim = c(min(c(Div_swap[3, ], Div_EffectSize[3, ])) - 5, 
              max(c(Div_swap[3, ], Div_EffectSize[3, ])) + 5),
     ylab = 'Diversity', xlab = 'Community')
lines(Div_EffectSize[3, ], col = 'red')


# Shared legend
par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  "center",
  legend = labs,
  col = cols,
  lty = 1,
  lwd = 2,
  horiz = TRUE,
  bty = "n"
)

# Reset layout
layout(1)
