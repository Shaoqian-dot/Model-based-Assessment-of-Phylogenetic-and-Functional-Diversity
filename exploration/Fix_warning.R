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
family_type <- c("poisson")
p <- c(80)
q <- 5
r <- 64
NOPS <- 16
alpha <- log(9/4)
beta <- log(2)
quantile <- 0.35
Corr <- 1
Eigen <- 1

###########################################################
## V_J and D_J
###########################################################
DM_phy_func <- getPhyloMatrix(tree = tree, m = p)
VC_phy_func <- 1 - DM_phy_func / max(DM_phy_func)
I_p <- diag(p)
J <- I_p - 1 / p * matrix(1, p, p)
S_J <- t(J) %*% VC_phy_func %*% J
eig <- spectral_decomp(VC_phy_func = S_J)
V_J <- eig$P
D_J <- eig$D

val_num.eig <- get_num.eig(
  Methods = 'Method3',
  Methods_para = NA,
  D = D_J
)
dat <- getData(
  DM_phy_func = DM_phy_func,
  NOPS = NOPS,
  q = q,
  alpha = alpha,
  beta = beta,
  r = r,
  quantile = quantile,
  Corr = Corr,
  P = V_J[, 1:(p - 1)],
  Eigen = Eigen,
  Distribution = family_type,
  p = p
)

dat_long <- dat$yX
dat_wide <- dat$abundance

###########################################################
## Modelling
###########################################################
fit_glmm_rr <- try(
  fit_model(
    type = "p_mid",
    matrix_type = "P_J",
    rr = TRUE,
    yX = dat_long,
    P_J = V_J,
    val_num.eig = val_num.eig,
    Distribution = family_type,
    null = FALSE,
    p = p
  ),
  silent = TRUE
)
fit_glmm_rr_fix <- try(
     fit_model(
         type = "p_mid",
         matrix_type = "P_J",
         rr = TRUE,
         yX = dat_long,
         P_J = V_J,
         val_num.eig = val_num.eig,
         Distribution = family_type,
         null = FALSE,
         p = p
       ),
     silent = TRUE
   )

cond <- rbind(
  a = fixef(fit_glmm_rr)$cond,
  b = fixef(fit_glmm_rr_fix)$cond
)
save()