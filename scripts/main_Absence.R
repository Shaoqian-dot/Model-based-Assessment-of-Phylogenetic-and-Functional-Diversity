# The threshold for rare species at sample size 64: 1/64 = 0.015625

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
library(stringr)
source("R/run_simulation_global.R")
source("R/simulation_core_global.R")
source("R/get_num.eig.R") # The default method is 'Method3'
source("R/fit_model_current2.R")
source("R/getData_current.R")
source("R/U_current.R")
source("R/getChisquare_fixed.R")
source("R/Test_LRT.R")
source("R/getChisquare.R")
source("R/getRaosQ.R")
source("R/RaoQ.R")
source("R/randomization.R")
source("R/rand.RaoQ.fun.R")
source('R/MeanAbund.R')
source('R/save_checkpoint.R')
# source('R/get_lrt_row.R')
# source('R/get_wald_row.R')

############################################################
### Settings
############################################################

set.seed(123)

############################################################
## Load parameters estimated from the PCV data set: tree, True_Para, Lambda
############################################################

load(file.path('data', "Para_PCV.RData"))

###########################################################
## Main simulation runner
###########################################################
axis_method <- 'Method3'
families <- c("binomial")
globalTest <- FALSE
test <- c("Both") # Or LRT
Swap = TRUE
p <- c(10)
q <- 5
r <- c(8)
alpha <- log(9/4)
beta <- log(2)
quantile <- 0.35
Eigen <- 3
NOS <- 2

# Other parameters from the previous simulations
signal <- NA
c_val <- NA
sigma2 <- NA
Corr <- 1

nsim = 1 # The number of simulations

seed <- data.frame(
  seed = seq(
    from = 1000,
    by = 1,
    length.out = nsim
  ))
res_all <- seed |>
  pmap_dfr(function(seed){
    
    run_simulation(
      family_type = families,
      globalTest = globalTest,
      test = test,
      Swap = Swap,
      p = p,
      q = q,
      r = r,
      seed = seed,
      tree = tree,
      alpha = alpha,
      beta = beta,
      quantile = quantile,
      Eigen = Eigen,
      axis_method = axis_method ,
      Absence = TRUE,
      FactorLoading = Lambda,
      TruePara = TruePara
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

