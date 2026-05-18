rm(list = ls())
############################################################
## Simulation study
## Compare powers of:
## 1. glmmTMB without rr
## 2. glmmTMB with rr
## 3. Mlm()
## 4. lm()
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

############################################################
### Settings
############################################################

set.seed(123)

m_vec  <- c(2, 4, 8, 16)
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

############################################################
### Function: Data generation
############################################################

generate_data <- function(m,
                          r,
                          signal = c("v1", "vm"),
                          c_val = 0.6,
                          sigma2 = 1){
  
  signal <- match.arg(signal)
  
  ##########################################################
  ## Phylogenetic similarity matrix and its eigenvectors
  ##########################################################
  
  Phy_DM <- getPhyloMatrix(tree, m)
  
  Phy_SM <- max(Phy_DM) - Phy_DM
  
  spec <- spectral_decomp(VC_phy_func = Phy_SM)
  
  V <- spec$P
  
  vec1 <- V[, 1]
  vecm <- V[, m]
  
  ##########################################################
  ## Factor loadings
  ##########################################################
  
  Lambda <- cbind(
    V[, 1],
    V[, 2]
  )
  
  ##########################################################
  ## Mean structure
  ##########################################################
  
  mu1 <- rep(0, m)
  
  if(signal == "v1"){
    mu2 <- c_val * vec1
  }
  
  if(signal == "vm"){
    mu2 <- c_val * vecm
  }
  
  ##########################################################
  ## Residual latent structure
  ##########################################################
  
  ## z_il ~ N(0, I_2)
  
  Z1 <- matrix(rnorm(r * 2), r, 2)
  Z2 <- matrix(rnorm(r * 2), r, 2)
  
  U1 <- Z1 %*% t(Lambda)
  U2 <- Z2 %*% t(Lambda)
  
  ##########################################################
  ## Generate observations
  ##########################################################
  
  Y1 <- matrix(0, r, m)
  Y2 <- matrix(0, r, m)
  
  for(j in 1:m){
    
    Y1[, j] <- mu1[j] + U1[, j] + rnorm(r, 0, sqrt(sigma2))
    
    Y2[, j] <- mu2[j] + U2[, j] + rnorm(r, 0, sqrt(sigma2))
  }
  
  ##########################################################
  ## Long-format data for glmmTMB
  ##########################################################
  
  Y <- rbind(Y1, Y2)
  
  dat_long <- data.frame(
    id = rep(1:(2 * r), each = m),
    com = rep(rep(c(1, 2), each = r), each = m),
    sp = factor(rep(1:m, times = 2 * r)),
    y = c(t(Y))
  )
  
  ##########################################################
  ## Wide-format data for lm / Mlm
  ##########################################################
  
  dat_wide <- data.frame(
    com = factor(rep(c(1, 2), each = r)),
    Y
  )
  
  colnames(dat_wide)[-1] <- paste0("sp", 1:m)
  
  return(list(
    long = dat_long,
    wide = dat_wide
  ))
}

############################################################
### Function: Fit models
############################################################

fit_models <- function(dat, m){
  
  dat_long <- dat$long
  dat_wide <- dat$wide
  
  ##########################################################
  ## glmmTMB without rr
  ##########################################################
  
  fit_glmm_no_rr <- try(
    
    glmmTMB(
      y ~ com * sp,
      data = dat_long,
      family = gaussian()
    ),
    
    silent = TRUE
  )
  
  fit_glmm_no_rr_null <- try(
    
    glmmTMB(
      y ~ sp,
      data = dat_long,
      family = gaussian()
    ),
    
    silent = TRUE
  )
  
  ##########################################################
  ## glmmTMB with rr
  ##########################################################
  
  fit_glmm_rr <- try(
    
    glmmTMB(
      y ~ com * sp +
        rr(sp + 0 | id, d = 2),
      data = dat_long,
      family = gaussian()
    ),
    
    silent = TRUE
  )
  
  fit_glmm_rr_null <- try(
    
    glmmTMB(
      y ~ sp +
        rr(sp + 0 | id, d = 2),
      data = dat_long,
      family = gaussian()
    ),
    
    silent = TRUE
  )
  
  ##########################################################
  ## mlm object
  ##########################################################
  
  Ymat <- as.matrix(dat_wide[, -1])
  
  fit_mlm <- try(
    
    lm(
      Ymat ~ com,
      data = dat_wide
    ),
    
    silent = TRUE
  )
  
  ##########################################################
  ## lm()
  ##########################################################
  
  pvals_lm <- rep(NA, m)
  
  for(j in 1:m){
    
    tmp_fit <- lm(
      dat_wide[[j + 1]] ~ dat_wide$com
    )
    
    pvals_lm[j] <- anova(tmp_fit)[1, "Pr(>F)"]
  }
  
  ##########################################################
  ## Likelihood ratio tests
  ##########################################################
  
  get_lrt_p <- function(fit1, fit0){
    
    if(inherits(fit1, "try-error") ||
       inherits(fit0, "try-error")){
      return(NA)
    }
    
    out <- try(
      anova(fit0, fit1),
      silent = TRUE
    )
    
    if(inherits(out, "try-error")){
      return(NA)
    }
    
    return(out$`Pr(>Chisq)`[2])
  }
  
  p_glmm_no_rr <- get_lrt_p(
    fit_glmm_no_rr,
    fit_glmm_no_rr_null
  )
  
  p_glmm_rr <- get_lrt_p(
    fit_glmm_rr,
    fit_glmm_rr_null
  )
  
  ##########################################################
  ## Mlm p-value
  ##########################################################
  
  p_mlm <- NA
  
  if(!inherits(fit_mlm, "try-error")){
    
    tmp <- try(
      summary(manova(fit_mlm)),
      silent = TRUE
    )
    
    if(!inherits(tmp, "try-error")){
      p_mlm <- tmp$stats['com', 'Pr(>F)']
    }
  }
  
  ##########################################################
  ## Combine lm p-values
  ##########################################################
  
  ## Fisher method
  
  p_lm <- pchisq(
    -2 * sum(log(pvals_lm)),
    df = 2 * m,
    lower.tail = FALSE
  )
  
  ##########################################################
  ## Return
  ##########################################################
  
  return(data.frame(
    glmm_no_rr = p_glmm_no_rr,
    glmm_rr    = p_glmm_rr,
    mlm         = p_mlm,
    lm          = p_lm
  ))
}

############################################################
### Simulation
############################################################

run_simulation <- function(signal_type){
  
  res_all <- list()
  
  counter <- 1
  
  for(m in m_vec){
    
    cat("Running: signal =", signal_type,
        ", m =", m, "\n")
    
    
    ########################################################
    ## Replications
    ########################################################
    
    tmp_res <- map_dfr(
      1:nsim,
      function(sim){
        
        dat <- generate_data(
          m = m,
          r = r,
          signal = signal_type,
          c_val = c_val,
          sigma2 = sigma2
        )
        
        fit_models(dat, m)
      }
    )
    
    ########################################################
    ## Power
    ########################################################
    
    power_res <- colMeans(
      tmp_res < 0.05,
      na.rm = TRUE
    )
    
    res_all[[counter]] <- data.frame(
      signal = signal_type,
      m = m,
      model = names(power_res),
      power = as.numeric(power_res)
    )
    
    counter <- counter + 1
  }
  
  bind_rows(res_all)
}

############################################################
### Run simulation
############################################################

res_v1 <- run_simulation("v1")

res_vm <- run_simulation("vm")

res_all <- bind_rows(
  res_v1,
  res_vm
)

saveRDS(
  res_all,
  file = file.path("results", "res_all.rds")
)

############################################################
### Power plot
############################################################

ggplot(
  res_all,
  aes(x = m,
      y = power,
      color = model)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ signal) +
  scale_x_continuous(
    breaks = m_vec
  ) +
  ylim(0, 1) +
  theme_bw() +
  labs(
    x = "Number of species (m)",
    y = "Power",
    color = "Model"
  )
