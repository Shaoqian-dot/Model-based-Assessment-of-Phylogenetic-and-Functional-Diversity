# rm(list = ls())
# 
# ############################################################
# ## Simulation runner for Katana HPC
# ##
# ## This script:
# ## 1. Reads command line arguments
# ## 2. Loads simulation functions
# ## 3. Runs one simulation job
# ## 4. Saves the output as an .rds file
# ##
# ## Intended usage:
# ##
# ## Rscript scripts/run_simulation.R \
# ##   family signal r c_val sigma2 globalTest p seed
# ##
# ############################################################
# 
# ############################################################
# ## Print runtime information
# ############################################################
# 
# cat("========================================\n")
# cat("Starting simulation job\n")
# cat("Working directory:", getwd(), "\n")
# 
# library(here)
# cat("Project root (here):", here(), "\n")
# cat("========================================\n\n")
# 
# ############################################################
# ## Load required packages
# ############################################################
# 
# library(MASS)
# library(glmmTMB)
# library(vegan)
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# library(purrr)
# library(ape)
# 
# 
# ############################################################
# ## Source custom functions
# ############################################################
# 
# source(here("R", "getPhyloMatrix.R"))
# source(here("R", "spectral_decomp.R"))
# source(here("R", "run_simulation_global.R"))
# source(here("R", "simulation_core_global.R"))
# source(here("R", "get_num.eig.R")) # The default method is 'Method3'
# source(here("R", "fit_model_current2.R"))
# source(here("R", "getData_current.R"))
# source(here("R", "GroupingPairs.R"))
# source(here("R", "mean_matrix_current.R"))
# source(here("R", "sample_pairs.R"))
# source(here("R", "U_current.R"))
# source(here("R", "getChisquare_fixed.R"))
# source(here("R", "getData_nonSwap.R"))
# source(here("R", "Test_LRT.R"))
# source(here("R", "getChisquare.R"))
# source(here("R", "getRaosQ.R"))
# source(here("R", "RaoQ.R"))
# source(here("R", "randomization.R"))
# source(here("R", "rand.RaoQ.fun.R"))
# 
# ############################################################
# ## Read simulation parameters
# ############################################################
# 
# ## ---------------------------------------------------------
# ## Local debugging mode
# ## Mimic command-line arguments by reading the first row
# ## of config/params.csv. This block can be used when
# ## testing the script interactively in RStudio or locally.
# ## ---------------------------------------------------------
# 
# # params <- read.csv("config/params.csv", stringsAsFactors = FALSE)
# # 
# # # Convert the first parameter set into a character vector
# # # matching the format returned by commandArgs()
# # args <- as.character(unlist(params[1, ]))
# # 
# # print(args)
# 
# ## ---------------------------------------------------------
# ## Katana / PBS execution mode
# ## Read arguments supplied from the job submission script.
# ## This overwrites the local debugging arguments above.
# ## ---------------------------------------------------------
# 
# args <- commandArgs(trailingOnly = TRUE)
# 
# expected_n <- 1
# 
# # Check that the expected number of arguments has been
# # supplied before proceeding with the simulation.
# if (length(args) != expected_n) {
#   stop(
#     paste0(
#       "Expected ", expected_n, " command line arguments.\n",
#       "Received: ", length(args)
#     )
#   )
# }
# 
# ############################################################
# ## Parse arguments
# ############################################################
# 
# family      <- c('binomial', 'poisson')
# signal      <- NA
# 
# p           <- C(5, 10, 20, 40, 80)
# r           <- c(4, 8, 16, 32, 64)
# c_val       <- NA
# sigma2      <- NA
# 
# globalTest  <- FALSE
# test        <- 'Both'
# Swap        <- 'TRUE'  
# alpha       <- log(9/4)
# beta        <- log(2)
# quantile    <- 0.35
# 
# Corr        <- 1
# Eigen       <- c(1, 2)
# Method      <- axis_method <- 'Both'
# 
# seed        <- as.numeric(args[1])
# 
# q <- 5
# nsim <- 1000
# 
# 
# ############################################################
# ## Set random seed
# ############################################################
# 
# set.seed(seed)
# 
# ############################################################
# ## Print simulation settings
# ############################################################
# 
# cat("Simulation settings:\n")
# cat("----------------------------------------\n")
# cat("family      :", family, "\n")
# cat("signal      :", signal, "\n")
# cat("r           :", r, "\n")
# cat("c_val       :", c_val, "\n")
# cat("sigma2      :", sigma2, "\n")
# cat("globalTest  :", globalTest, "\n")
# cat("p           :", p, "\n")
# cat("seed        :", seed, "\n")
# cat("test        :", test, "\n")
# cat("Swap        :", Swap, "\n")
# cat("alpha       :", alpha, "\n")
# cat("beta        :", beta, "\n")
# cat("quantile    :", quantile, "\n")
# cat("Corr        :", Corr, "\n")
# cat("Eigen       :", Eigen, "\n")
# cat("Method      :", Method, "\n")
# cat("----------------------------------------\n\n")
# 
# ############################################################
# ## Load phylogenetic tree
# ############################################################
# 
# tree_file <- here("data", "example.tre")
# 
# if(!file.exists(tree_file)){
#   stop("Tree file not found: ", tree_file)
# }
# 
# tree <- ape::read.tree(tree_file)
# 
# ############################################################
# ## Run simulation
# ############################################################
# 
# cat("Running simulation...\n\n")
# 
# res <- run_simulation(
#   family_type = family,
#   signal_type = signal,
#   r = r,
#   c_val = c_val,
#   sigma2 = sigma2,
#   globalTest = globalTest,
#   p = p,
#   q = q,
#   seed = seed,
#   tree = tree,
#   test = test,
#   Swap = Swap,
#   alpha = alpha,
#   beta = beta,
#   quantile = quantile,
#   Corr = Corr,
#   Eigen = Eigen,
#   Method = Method
# )
# 
# ############################################################
# ## Create output directory
# ############################################################
# 
# run_name <- "Diversity_Poisson_p5_Both_Method3_S2_RE(p-k)_trail"
# # S1 : 1 - DM_phy_func / (max(DM_phy_func) + 1)
# # S2 : 1 - DM_phy_func /max(DM_phy_func)
# # RE : the dimension of multivariate random effect
# dir.create(
#   here("output", run_name),
#   recursive = TRUE,
#   showWarnings = FALSE
# )
# 
# ############################################################
# ## Define output filename
# ############################################################
# 
# outfile <- here(
#   "output",
#   run_name,
#   paste0(
#     "res_seed",
#     seed,
#     ".rds"
#   )
# )
# 
# ############################################################
# ## Save results
# ############################################################
# 
# saveRDS(res, outfile)
# 
# ############################################################
# ## Finish message
# ############################################################
# 
# cat("\nSimulation completed successfully.\n")
# cat("Output saved to:\n")
# cat(outfile, "\n")
# cat("========================================\n")
# 
# 



rm(list = ls())

############################################################
## Simulation runner for Katana HPC
##
## This script
##   1. Reads the job seed from the command line.
##   2. Loads all simulation functions.
##   3. Runs one simulation job.
##   4. Saves the results as an .rds file.
##
## Usage:
##
## Rscript scripts/run_simulation.R <seed>
##
############################################################

############################################################
## Runtime information
############################################################

cat("========================================\n")
cat("Starting simulation job\n")
cat("Working directory:", getwd(), "\n")

library(here)

cat("Project root:", here(), "\n")
cat("========================================\n\n")

############################################################
## Packages
############################################################

library(MASS)
library(glmmTMB)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ape)

############################################################
## Source functions
############################################################

source(here("R", "getPhyloMatrix.R"))
source(here("R", "spectral_decomp.R"))
source(here("R", "run_simulation_global.R"))
source(here("R", "simulation_core_global.R"))
source(here("R", "get_num.eig.R"))
source(here("R", "fit_model_current2.R"))
source(here("R", "getData_current.R"))
source(here("R", "GroupingPairs.R"))
source(here("R", "mean_matrix_current.R"))
source(here("R", "sample_pairs.R"))
source(here("R", "U_current.R"))
source(here("R", "getChisquare_fixed.R"))
source(here("R", "getData_nonSwap.R"))
source(here("R", "Test_LRT.R"))
source(here("R", "getChisquare.R"))
source(here("R", "getRaosQ.R"))
source(here("R", "RaoQ.R"))
source(here("R", "randomization.R"))
source(here("R", "rand.RaoQ.fun.R"))

############################################################
## Read command line arguments
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript run_simulation.R <seed>")
}

seed <- as.numeric(args[1])

############################################################
## Simulation settings
############################################################

family <- c("poisson")

p <- c(80)

r <- c(64)

Eigen <- c(1)


q <- 5

signal <- NA

c_val <- NA
sigma2 <- NA

globalTest <- FALSE

test <- "Both"

Swap <- TRUE

alpha <- log(9 / 4)
beta <- log(2)
quantile <- 0.35

Corr <- 1

axis_method <- "Both"

############################################################
## Set random seed
############################################################

set.seed(seed)

############################################################
## Print settings
############################################################

cat("Simulation settings\n")
cat("----------------------------------------\n")
cat("Seed         :", seed, "\n")
cat("Families     :", paste(family, collapse = ", "), "\n")
cat("p            :", paste(p, collapse = ", "), "\n")
cat("r            :", paste(r, collapse = ", "), "\n")
cat("Eigen        :", paste(Eigen, collapse = ", "), "\n")
cat("Method       :", Method, "\n")
cat("Global test  :", globalTest, "\n")
cat("Test         :", test, "\n")
cat("Swap         :", Swap, "\n")
cat("----------------------------------------\n\n")

############################################################
## Load tree
############################################################

tree_file <- here("data", "example.tre")

if (!file.exists(tree_file)) {
  stop("Cannot find tree file: ", tree_file)
}

tree <- ape::read.tree(tree_file)

############################################################
## Run simulation
############################################################

cat("Running simulations...\n\n")

res <- run_simulation(
  signal_type = signal,
  family_type = family,
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
  axis_method  = axis_method 
)

############################################################
## Save results
############################################################

run_name <- "Diversity_Pois_p_r_TestBoth_MethodBoth_S2_RE(p-k)_trial"

out_dir <- here("output", run_name)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

outfile <- file.path(
  out_dir,
  paste0("res_seed", seed, ".rds")
)

saveRDS(res, outfile)

############################################################
## Finish
############################################################

cat("\nSimulation completed successfully.\n")
cat("Results saved to:\n")
cat(outfile, "\n")
cat("========================================\n")