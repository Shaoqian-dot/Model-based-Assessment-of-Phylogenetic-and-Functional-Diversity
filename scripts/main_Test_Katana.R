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

family <- c("binomial")

p <- c(5, 10, 20, 40, 80)

r <- c(4, 8, 16, 32, 64)

Eigen <- c(2)


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

run_name <- "PhyD_Bin_p510204080_r48163264_TestBoth_MethodBoth_Eigen2"

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