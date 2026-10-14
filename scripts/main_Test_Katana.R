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
library(stringr)

############################################################
## Source functions
############################################################

source(here("R", "getPhyloMatrix.R"))
source(here("R", "spectral_decomp.R"))
source(here("R", "run_simulation_global.R"))#
source(here("R", "simulation_core_global.R"))#
source(here("R", "get_num.eig.R"))#
source(here("R", "fit_model_current2.R"))#
source(here("R", "getData_current.R"))#
source(here("R", "GroupingPairs.R"))
source(here("R", "mean_matrix_current.R"))
source(here("R", "sample_pairs.R"))
source(here("R", "U_current.R"))#
source(here("R", "getChisquare_fixed.R"))#
source(here("R", "getData_nonSwap.R"))
source(here("R", "Test_LRT.R"))#
source(here("R", "getChisquare.R"))#
source(here("R", "getRaosQ.R"))#
source(here("R", "RaoQ.R"))#
source(here("R", "randomization.R"))#
source(here("R", "rand.RaoQ.fun.R"))#
source(here("R", "MeanAbund.R"))#
source(here("R", "save_checkpoint.R"))#

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

family <- c("binomial")

p <- c(5, 10, 20, 40, 80)

r <- c(4, 8, 16, 32, 64)

Eigen <- c(3)


q <- 5

# signal <- NA
# 
# c_val <- NA
# sigma2 <- NA

#Corr <- 1

globalTest <- FALSE

test <- "Both"

Swap <- TRUE

alpha <- log(9 / 4)
beta <- log(4)
quantile <- 0.35
axis_method <- "Method3"
Absence <- TRUE

############################################################
## Set random seed
############################################################

set.seed(seed)

############################################################
## Print settings
############################################################

cat("Simulation settings\n")
cat("----------------------------------------\n")
cat("Seed              :", seed, "\n")
cat("Families          :", paste(family, collapse = ", "), "\n")
cat("p                 :", paste(p, collapse = ", "), "\n")
cat("r                 :", paste(r, collapse = ", "), "\n")
cat("Eigen             :", paste(Eigen, collapse = ", "), "\n")
cat("axis_method       :", axis_method, "\n")
cat("Global test       :", globalTest, "\n")
cat("Test              :", test, "\n")
cat("Swap              :", Swap, "\n")
cat("----------------------------------------\n\n")

############################################################
## Load parameters estimated from the PCV data set: tree, True_Para, Lambda
############################################################

load(file.path('data', "Para_PCV.RData"))

############################################################
## Run simulation
############################################################

cat("Running simulations...\n\n")

res <- run_simulation(
  family_type = family,
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
  axis_method  = axis_method,
  Absence = Absence,
  FactorLoading = Lambda,
  TruePara = TruePara
)

############################################################
## Save results
############################################################

run_name <- paste0('Phy_', 'fam_', family,'_test_', test,
                   '_axis_method_', axis_method, '_Eigen_', Eigen, '_Absence_', Absence)

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