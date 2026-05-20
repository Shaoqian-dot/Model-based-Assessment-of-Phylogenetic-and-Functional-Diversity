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
source("R/get_num.eig.R")
source("R/fit_model_current2.R")
############################################################
### Settings
############################################################

set.seed(123)

p_vec  <- c(2, 4, 8, 16)
r      <- 80
nsim   <- 100
sigma2 <- 1
c_val  <- 0.6

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

families <- c("poisson")
signals <- c("v1", "vp")

res_all <- expand.grid(family = families,
                       signal = signals,
                       stringsAsFactors = FALSE) |>
  pmap_dfr(function(family, signal){
    
    run_simulation(signal = signal,
                   family = family, 
                   globalTest = FALSE)
  })

saveRDS(
  res_all,
  file = file.path("results", "res_all_Poisson_nonGlobal.rds")
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

