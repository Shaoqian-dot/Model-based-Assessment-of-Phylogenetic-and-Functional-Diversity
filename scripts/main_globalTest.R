rm(list = ls())
############################################################
## Simulation study
## Compare powers of:
## 1. glmmTMB without rr
## 2. glmmTMB with rr
############################################################

library(MASS)
library(glmmTMB)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ape)
source("R/getPhyloMatrix.R")
source("R/spectral_decomp.R")
source("R/run_simulation_global.R")
source("R/simulation_core_global.R")
source("R/get_num.eig.R") # The default method is 'Method3'
source("R/fit_model_current2.R")
source("R/getData_current.R")
source("R/GroupingPairs.R")
source("R/mean_matrix_current.R")
source("R/sample_pairs.R")
source("R/U_current.R")
source("R/getChisquare_fixed.R")
source("R/getData_nonSwap.R")
source("R/Test_LRT.R")
source("R/getChisquare.R")
source("R/getRaosQ.R")
source("R/RaoQ.R")
source("R/randomization.R")
source("R/rand.RaoQ.fun.R")
#source("R/Test_Wald.R")
############################################################
### Settings
############################################################

set.seed(123)

############################################################
### Load phylogenetic tree
############################################################

data_dir <- "data"

tree_file <- file.path(data_dir, "example.tre")

if(!file.exists(tree_file)){
  stop("Tree file not found: ", tree_file)
}

tree <- ape::read.tree(tree_file)

###########################################################
## Main simulation runner
###########################################################

axis_method <- 'Method3'
signal <- NA
families <- c("poisson")
globalTest <- FALSE
test <- c("Both") # Or LRT
Swap = TRUE
p <- c(5)
q <- 5
r <- c(64)
c_val <- NA
sigma2 <- NA
alpha <- log(9/4)
beta <- log(2)
quantile <- 0.35
Corr <- 1
Eigen <- 1


nsim = 1 # The number of simulations
step <- length(p) * length(r) * length(families) * length(Eigen)

seed <- data.frame(
  seed = seq(
  from = 1000,
  by = step,
  length.out = nsim
))
res_all <- seed |>
  pmap_dfr(function(seed){
    
    run_simulation(
      signal_type = signal,
      family_type = families,
      globalTest = globalTest,
      test = test,
      Swap = Swap,
      p = p,
      q = q,
      r = r,
      c_val = c_val,
      sigma2 = sigma2,
      seed = seed,
      tree = tree,
      alpha = alpha,
      beta = beta,
      quantile = quantile,
      Corr = Corr,
      Eigen = Eigen,
      axis_method = axis_method 
    )
  })

res_mean <- res_all %>%
  group_by(family,
           signal,
           p,
           r,
           c_val,
           model,
           globalTest,
           test,
           Eigen,
           com,
           axis_method,
           method) %>%
  summarise(
    power = mean(p_values < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(
  res_all,
  file = file.path("results", "res_all_Pois_NegBino_Global_c0.6.rds")
)

############################################################
### Power plot
############################################################

ggplot(
  res_all,
  aes(x = p,
      y = power,
      color = model)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ signal) +
  scale_x_continuous(
    breaks = p_vec
  ) +
  ylim(0, 1) +
  theme_bw() +
  labs(
    x = "Number of species (p)",
    y = "Power",
    color = "Model"
  )

