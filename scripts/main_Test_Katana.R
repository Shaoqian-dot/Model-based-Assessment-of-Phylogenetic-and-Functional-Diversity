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

source(here("R", "getPhyloMatrix.R"))
source(here("R", "spectral_decomp.R"))
source(here("R", "run_simulation_global.R"))
source(here("R", "simulation_core_global.R"))
source(here("R", "get_num.eig.R")) # The default method is 'Method3'
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

############################################################
## Read simulation parameters
############################################################

## ---------------------------------------------------------
## Local debugging mode
## Mimic command-line arguments by reading the first row
## of config/params.csv. This block can be used when
## testing the script interactively in RStudio or locally.
## ---------------------------------------------------------

# params <- read.csv("config/params.csv", stringsAsFactors = FALSE)
# 
# # Convert the first parameter set into a character vector
# # matching the format returned by commandArgs()
# args <- as.character(unlist(params[1, ]))
# 
# print(args)

## ---------------------------------------------------------
## Katana / PBS execution mode
## Read arguments supplied from the job submission script.
## This overwrites the local debugging arguments above.
## ---------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

expected_n <- 18

# Check that the expected number of arguments has been
# supplied before proceeding with the simulation.
if (length(args) != expected_n) {
  stop(
    paste0(
      "Expected ", expected_n, " command line arguments.\n",
      "Received: ", length(args)
    )
  )
}

############################################################
## Parse arguments
############################################################

family      <- args[1]
signal      <- args[2]

r           <- as.numeric(args[3])
c_val       <- as.numeric(args[4])
sigma2      <- as.numeric(args[5])

globalTest  <- as.logical(args[6])
test        <- args[7]  
Swap        <- as.logical(args[8])  
alpha       <- as.numeric(args[9])
beta        <- as.numeric(args[10])
quantile    <- as.numeric(args[11]) 

Corr        <- as.numeric(args[12])
Eigen       <- as.numeric(args[13])

p           <- as.numeric(args[14])
nsim        <- as.numeric(args[15])
seed        <- as.numeric(args[16])
Method      <- args[17]
Methods_para<- as.numeric(args[18])
q <- 5





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
cat("test        :", test, "\n")
cat("Swap        :", Swap, "\n")
cat("alpha       :", alpha, "\n")
cat("beta        :", beta, "\n")
cat("quantile    :", quantile, "\n")
cat("Corr        :", Corr, "\n")
cat("Eigen       :", Eigen, "\n")
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
  q = q,
  nsim = nsim,
  seed = seed,
  tree = tree,
  test = test,
  Swap = Swap,
  alpha = alpha,
  beta = beta,
  quantile = quantile,
  Corr = Corr,
  Eigen = Eigen,
  Method = Method,
  Methods_para = Methods_para
)

############################################################
## Create output directory
############################################################

run_name <- "Diversity_Poisson_p5_Wald_Method1_2"

dir.create(
  here("output", run_name),
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
## Define output filename
############################################################

outfile <- here(
  "output",
  run_name,
  paste0(
    "res_seed",
    seed,
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