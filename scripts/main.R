rm(list = ls())
options(error = recover)  # Enable interactive debugging on error

#browser()
#start_time <- Sys.time()

# This script performs simulation using the first eigenvector of a distance matrix as contrast

# -------------------- PACKAGE SETUP --------------------
# install.packages("scales")
# install.packages('vegan')
# install.packages("glmmTMB")
# install.packages('multcomp')
# install.packages('ape')
# install.packages('picante')
# install.packages('plotly')
# install.packages('patchwork')
# install.packages('FD')
# install.packages('progress')
# install.packages('tidyverse')
# install.packages('DHARMa')

library(glmmTMB)
library(multcomp)
library(ape)
library(picante)        # For tipShuffle()
library(MASS)
library(scales)
library(ggplot2)
library(tidyverse)
library(plotly)
library(htmlwidgets)
library(tidyr)
library(dplyr)
library(patchwork)      # Plot composition
library(FD)             # Functional diversity (gowdis)
library(tibble)
library(purrr)
library(progress)
library(Matrix)
library(DHARMa)         # Residual diagnostics
library(ecostats)
library(vegan)
library(emmeans) 
################################################################################
######### MAIN SIMULATION SCRIPT (EXECUTED VIA PBS JOB ARRAY ON KATANA) #########
################################################################################

################################################################################
# -------------------- PATH SETUP --------------------
setwd('C:/Users/huang/OneDrive - UNSW/PhD Project/Model contrasts/Simulation/Simulation with wide settings/Code refactoring/Model-based assessment of PD and FD/Katana/Assess_Diversity_Corr_FL_P_current_negBino_LRT')
this_dir <- getwd()
data_dir <- file.path(this_dir, "data")

# -------------------- LOAD CUSTOM FUNCTIONS --------------------
source('functions/PD_FDSimulation_Correlation_4.R')
source('functions/mean_matrix_current.R')
source('functions/rand.RaoQ.fun.R')
source('functions/randomization.R')
source('functions/RaoQ.R')
source('functions/GroupingPairs.R')
source('functions/getFuncMatrix.R')
source('functions/getPhyloMatrix.R')
source('functions/spectral_decomp.R')
source('functions/getData_current.R')
source('functions/getRaosQ.R')
source('functions/get_num.eig.R')
source('functions/getChisquare_fixed.R')
source('functions/getChisquare.R')
source('functions/Ratio_NA_NaN_calculation.R')
source('functions/duplicate_random.R')
source('functions/U_current.R')
source('functions/nameContrast.R')
source('functions/sample_pairs.R')
source('functions/fit_model_current2.R')
source('functions/trun.R')
source('functions/Power_Calculation.R')
################################################################################
# -------------------- JOB ARRAY CONFIGURATION --------------------
tab <- read.csv("job_array_trail.csv")
n <- length(tab$job_number)
for (i in seq_along(tab$job_number)) {
  i <- 181
  # Check the progress
  job <- tab$job_number[i]
  cat("Processing:", i, "/", n, "\r")
  
  # Retrieve job index from PBS environment
  #job = as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
  # -------------------- SCENARIO PARAMETERS --------------------
  m        = tab$varying_parameter1[tab$job_number == job]  # Number of species
  r        = tab$varying_parameter2[tab$job_number == job]  # Sample size per community
  base_seed= tab$varying_parameter3[tab$job_number == job]  # Base seed for reproducibility
  diversity= tab$varying_parameter4[tab$job_number == job]  # 1 = PD, 2 = FD
  quantile = tab$varying_parameter5[tab$job_number == job]  # Threshold for distance partitioning
  Corr     = tab$varying_parameter6[tab$job_number == job]  # Correlation level
  Eigen    = tab$varying_parameter7[tab$job_number == job]  # 1: P[c(1, 2)]; 2: P[, c(m - 1, m)]
  rep      = tab$varying_parameter8[tab$job_number == job]
  
  Distribution = tab$Distribution[tab$job_number == job]
  scen     = tab$scenario[tab$job_number == job]            # Scenario label
  dir      = tab$save_location[tab$job_number == job]       # Output directory
  
  
  # -------------------- INPUT DATA --------------------
  tree_file <- file.path(data_dir, "example.tre")
  if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)
  
  tree      <- read.tree(tree_file)
  trait     <- read.csv(file.path(data_dir, "trait data.csv"))
  my.sample <- read.csv(file.path(data_dir, 'community data.csv'))
  
  ################################################################################
  # -------------------- BASELINE MODEL --------------------
  # Used to parameterize simulation inputs
  # model_summary <- readRDS(file.path(data_dir, "model_summary.rds"))
  # names(model_summary)[2] <- "Lambda"
  
  # -------------------- GLOBAL PARAMETERS --------------------
  q        <- 5        # Number of communities
  alpha    <- c(log(9/4))   # Alpha diversity effect
  beta     <- c(log(2))     # Beta diversity effect
  n_Model  <- 2             # Number of fitted models
  
  # Eigenvector selection methods
  Methods       <- c('Method1', 'Method3')
  Methods_para  <- c(2, NA)
  Method_label <- c("EVs 1–2", "λ > mean")
  Model_label <- c("p_mid_J", "p_mid_J_rr")
  # -------------------- TRAIT PREPROCESSING --------------------
  attach(trait)
  
  # Log-transform and standardize selected traits
  trait$height_mean         <- as.numeric(scale(log(height_mean))) 
  trait$X1000.S             <- as.numeric(scale(log(X1000.S))) 
  trait$SLA_mean            <- as.numeric(scale(log(SLA_mean))) 
  trait$flowering.duration  <- as.numeric(scale(log(flowering.duration)))
  
  detach(trait)
  
  # Remove collinear variables
  cols_to_drop <- c('flowering_start.7.9', 'LF_E',
                    'seed.disp._ZO', 'pollination_SELF', 'reprod_veg')
  
  trait <- trait %>% dplyr :: select(-all_of(cols_to_drop))
  
  # Replace missing values with 0
  trait[, c(7:19)][is.na(trait[, c(7:19)])] <- 0
  
  # Convert categorical traits to ordered factors
  for (j in c(7:19)){
    trait[, j] <- ordered(trait[, j])
  }
  
  # Align species names
  rownames(trait) <- colnames(my.sample[,3:298])
  
  ################################################################################
  # -------------------- PARAMETER DERIVATION --------------------
  # Number of unique species used in FD
  m_less <- switch(as.character(m),
                   "5"  = 3,
                   "10" = 6,
                   "20" = 15,
                   "40" = 24,
                   "80" = 48,
                   "160" = 96,
                   stop("Unknown m"))
  
  # Number of swap operations
  NOPS <- switch(as.character(m),
                 "5"  = 1,
                 "10" = 2,
                 "20" = 4,
                 "40" = 8,
                 "80" = 16,
                 "160" = 32,
                 stop("Unknown m"))
  
  ################################################################################
  # -------------------- INITIALIZE STORAGE --------------------
  p_Model_Array_4dim <- array(NA, c(q-1, n_Model, length(Methods), rep))
  
  dimnames(p_Model_Array_4dim) <- list(
    row    = paste0("com", 2:q),      # Communities
    column = Model_label,
    face   = Method_label,            # Eigenvector methods
    rep    = paste0("rep", 1:rep)
  )
  
  # Rao's Q results
  p_RaosQ_Matrix      <- matrix(NA, (q - 1), rep)
  p_RaosQ_rand_Matrix <- matrix(NA, (q - 1), rep)
  
  dimnames(p_RaosQ_Matrix) <- list(
    row = paste0("com", 2:q),
    column = paste0("rep", 1:rep)
  )
  
  dimnames(p_RaosQ_rand_Matrix) <- list(
    row = paste0("com", 2:q),
    column = paste0("rep", 1:rep)
  )
  
  # Ratio of sigular contrast matrix
  NA_contrast_Array <- array(NA, c(length(Methods), n_Model, rep))
  dimnames(NA_contrast_Array) <- list(
    row    = Method_label,
    column = Model_label,
    face   = paste0("rep", 1:rep)
  )
  
  ################################################################################
  # -------------------- MAIN SIMULATION LOOP --------------------
  for (j in 1 : rep){
    set.seed(base_seed + j)  # Ensure reproducibility per replicate
    
    Results <- PD_FDSimulation(
      alpha = alpha, beta = beta,
      q = q, m = m, m_less = m_less,
      NOPS = NOPS,
      r = r, diversity = diversity,
      Methods = Methods,
      Methods_para = Methods_para,
      tree = tree, trait = trait,
      quantile = quantile,
      Corr = Corr,
      n_Model = n_Model,
      Eigen = Eigen,
      Distribution = Distribution
    )
    
    # Store outputs
    p_Model_Array_4dim[, , , j] <- Results$p_Model_Array
    p_RaosQ_Matrix[, j]         <- Results$p_RaosQ
    p_RaosQ_rand_Matrix[, j]    <- Results$p_RaosQ_rand
    NA_contrast_Array[ , , j]   <- Results$NA_contrast
  }
  
  ################################################################################
  # -------------------- NA / NaN DIAGNOSTICS --------------------
  Ratio_NA_NaN_Model      <- apply(p_Model_Array_4dim, c(1, 2, 3), Ratio_NA_NaN_calculation)
  Ratio_NA_NaN_RaosQ      <- apply(p_RaosQ_Matrix, 1, Ratio_NA_NaN_calculation)
  Ratio_NA_NaN_RaosQ_rand <- apply(p_RaosQ_rand_Matrix, 1, Ratio_NA_NaN_calculation)
  Ratio_NA_contrast       <- apply(NA_contrast_Array, c(1, 2), mean)
  Ratio_NA_NaN = list(
    Model = Ratio_NA_NaN_Model,
    Contrast = Ratio_NA_contrast,
    RaosQ = Ratio_NA_NaN_RaosQ,
    RaosQ_rand = Ratio_NA_NaN_RaosQ_rand
  )
  
  ################################################################################
  # -------------------- STORE RESULTS --------------------
  Summary_p_Matrix <- tibble(
    alpha = alpha,
    beta  = beta,
    r     = r,
    m     = m,
    p_Model_Array_4dim   = list(p_Model_Array_4dim),
    p_RaosQ_Matrix       = list(p_RaosQ_Matrix),
    p_RaosQ_rand_Matrix  = list(p_RaosQ_rand_Matrix),
    Ratio_NA_NaN_Matrix  = list(Ratio_NA_NaN)
  )
  
  # Add metadata for later aggregation
  res <- add_column(
    Summary_p_Matrix,
    scenario = scen,
    sim      = base_seed,
    job      = job,
    Corr     = Corr,
    Eigen    =Eigen,
    Distribution = Distribution
  )
  
  # Save result (one file per job)
  #save(list = "res", file = paste0(dir, "/res_", job, ".RDATA"))
}

################################################################################
# -------------------- FUNCTION DESCRIPTION --------------------
# PD_FDSimulation()
#
# Inputs:
# alpha        : strength of alpha diversity change
# beta         : strength of beta diversity change
# q            : number of communities
# m            : number of focal species
# m_less       : number of unique species (FD)
# NOPS         : number of swap operations
# r            : sample size per community
# Methods      : eigenvector selection methods
# Methods_para : parameters for each method
# trait        : species trait matrix
# tree         : phylogenetic tree
#
# Outputs:
# p_Model_Array    : p-values from model-based approach
# p_RaosQ          : p-values from Rao's Q
# p_RaosQ_rand     : p-values from randomized Rao's Q