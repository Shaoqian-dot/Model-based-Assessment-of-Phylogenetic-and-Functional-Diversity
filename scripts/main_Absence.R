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
Habitat <- 'OA' # The target habitat. 8 represents OA. The habitat order is CG, LF, MF, NE, OA, OP, SE, SF.
# Ratio of rare species: 10/115 =0.086957 

############################################################
### Load the model and choose community mean abundances
############################################################
# Load the model
temp_env <- new.env()

load("results/Real Application/my_environment.RData", envir = temp_env)

model <- temp_env$model_SiteEffect

rm(temp_env)

p_total <- nlevels(model$frame$sp)

# eta_T_hat
nd <- model$frame
nd$id[] <- NA
nd$site[] <- NA
eta_T_hat <- predict(model, 
                     newdata = nd,
                     type = "link",
                     re.form = NULL)

Link_data <- data.frame(
  Community = model$frame$com,
  Link = eta_T_hat
)
boxplot(Link ~ Community, data = Link_data) # ChatGPT suggests using 

eta_T_hat <- matrix(
  as.numeric(eta_T_hat),
  ncol = p_total,
  byrow = FALSE
)
eta_T_hat_unique <- unique(eta_T_hat)
rownames(eta_T_hat_unique) <- levels(model$frame$com)
  
TruePara <- plogis(eta_T_hat_unique[Habitat, ])
names(TruePara) <- levels(model$frame$sp)
# Factor Loading
rep <- model$obj$env$report()
Lambda <- rep$fact_load[[3]]
rownames(Lambda) <- levels(model$frame$sp)

############################################################
### Load phylogenetic tree and community data
############################################################

my.sample <- read.csv(
  "data/community data.csv",
  stringsAsFactors = FALSE)

tree <- read.tree('data/example.tre')

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

# Save for Katana running
output_folder <- "data"
save(
  TruePara,
  tree,
  Lambda,
  file = file.path(output_folder, "Para_PCV.RData")
) 
###########################################################
## Main simulation runner
###########################################################
axis_method <- 'Method3'
families <- c("binomial")
globalTest <- FALSE
test <- c("Both") # Or LRT
Swap = TRUE
p <- c(5)
q <- 5
r <- c(64)
alpha <- log(9/4)
beta <- log(2)
quantile <- 0.35
Eigen <- 3

# Other parameters from the previous simulations
signal <- NA
c_val <- NA
sigma2 <- NA
Corr <- 1

nsim = 2 # The number of simulations
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

