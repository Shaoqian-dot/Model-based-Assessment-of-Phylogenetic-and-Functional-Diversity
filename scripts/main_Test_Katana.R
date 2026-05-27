rm(list = ls())

cat("Working directory:", getwd(), "\n")
cat("here() root:", here(), "\n")

library(MASS)
library(glmmTMB)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(ape)
library(here)

############################################################
## Source functions
############################################################

source(here("R", "functions", "getPhyloMatrix.R"))
source(here("R", "functions", "spectral_decomp.R"))
source(here("R", "functions", "run_simulation_global.R"))
source(here("R", "functions", "simulation_core_global.R"))
source(here("R", "functions", "get_num.eig.R"))
source(here("R", "functions", "fit_model_current2.R"))

############################################################
## Command line arguments
############################################################

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 9){
  stop("Expected 9 command line arguments.")
}

family <- args[1]
signal <- args[2]
r <- as.numeric(args[3])
c_val <- as.numeric(args[4])
sigma2 <- as.numeric(args[5])
globalTest <- as.logical(args[6])
p <- as.numeric(args[7])
nsim   <- as.numeric(args[8])
seed   <- as.numeric(args[9])

set.seed(seed)

############################################################
## Load phylogenetic tree
############################################################

tree_file <- here("data", "example.tre")

if(!file.exists(tree_file)){
  stop("Tree file not found: ", tree_file)
}

tree <- ape::read.tree(tree_file)

############################################################
## Run simulation
############################################################

res <- run_simulation(
  family_type = family,
  signal_type = signal,
  r = r,
  c_val = c_val,
  sigma2 = sigma2,
  globalTest = globalTest,
  p = p,
  nsim = nsim,
  seed = seed,
  tree = tree
)

############################################################
## Save output
############################################################

outfile <- here(
  "output",
  paste0(
    "res_",
    family, "_",
    signal, "_",
    p, "_",
    c_val, "_",
    seed,
    ".rds"
  )
)

dir.create(here("output"), showWarnings = FALSE)
saveRDS(res, outfile)

cat("Finished:", outfile, "\n")

# params <- expand.grid(
#   signal = c("v1", "vp"),
#   family = c("poisson"),
#   rep = 1:100,
#   c_val = c(0.3, 0.6, 0.9)
# )
# params$seed <- 1000 + seq_len(nrow(params))




rm(list = ls())

############################################################
## Simulation runner for Katana HPC
##
## This script:
## 1. Reads command line arguments
## 2. Loads simulation functions
## 3. Runs one simulation job
## 4. Saves the output as an .rds file
##
## Intended usage:
##
## Rscript scripts/run_simulation.R \
##   family signal r c_val sigma2 globalTest p nsim seed
##
############################################################

############################################################
## Print runtime information
############################################################

cat("========================================\n")
cat("Starting simulation job\n")
cat("Working directory:", getwd(), "\n")

library(here)
cat("Project root (here):", here(), "\n")
cat("========================================\n\n")

############################################################
## Load required packages
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
## Source custom functions
############################################################

source(here("R", "functions", "getPhyloMatrix.R"))
source(here("R", "functions", "spectral_decomp.R"))
source(here("R", "functions", "run_simulation_global.R"))
source(here("R", "functions", "simulation_core_global.R"))
source(here("R", "functions", "get_num.eig.R"))
source(here("R", "functions", "fit_model_current2.R"))

############################################################
## Read command line arguments
############################################################

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 9){
  stop(
    paste0(
      "Expected 9 command line arguments.\n",
      "Received: ", length(args)
    )
  )
}

############################################################
## Parse arguments
############################################################

family     <- args[1]
signal     <- args[2]
r           <- as.numeric(args[3])
c_val       <- as.numeric(args[4])
sigma2      <- as.numeric(args[5])
globalTest <- as.logical(args[6])
p           <- as.numeric(args[7])
nsim        <- as.numeric(args[8])
seed        <- as.numeric(args[9])

############################################################
## Set random seed
############################################################

set.seed(seed)

############################################################
## Print simulation settings
############################################################

cat("Simulation settings:\n")
cat("----------------------------------------\n")
cat("family      :", family, "\n")
cat("signal      :", signal, "\n")
cat("r           :", r, "\n")
cat("c_val       :", c_val, "\n")
cat("sigma2      :", sigma2, "\n")
cat("globalTest  :", globalTest, "\n")
cat("p           :", p, "\n")
cat("nsim        :", nsim, "\n")
cat("seed        :", seed, "\n")
cat("----------------------------------------\n\n")

############################################################
## Load phylogenetic tree
############################################################

tree_file <- here("data", "example.tre")

if(!file.exists(tree_file)){
  stop("Tree file not found: ", tree_file)
}

tree <- ape::read.tree(tree_file)

############################################################
## Run simulation
############################################################

cat("Running simulation...\n\n")

res <- run_simulation(
  family_type = family,
  signal_type = signal,
  r = r,
  c_val = c_val,
  sigma2 = sigma2,
  globalTest = globalTest,
  p = p,
  nsim = nsim,
  seed = seed,
  tree = tree
)

############################################################
## Create output directory
############################################################

dir.create(
  here("output"),
  showWarnings = FALSE,
  recursive = TRUE
)

############################################################
## Define output filename
############################################################

outfile <- here(
  "output",
  paste0(
    "res_",
    family, "_",
    signal, "_p", p,
    "_c", c_val,
    "_seed", seed,
    ".rds"
  )
)

############################################################
## Save results
############################################################

saveRDS(res, outfile)

############################################################
## Finish message
############################################################

cat("\nSimulation completed successfully.\n")
cat("Output saved to:\n")
cat(outfile, "\n")
cat("========================================\n")