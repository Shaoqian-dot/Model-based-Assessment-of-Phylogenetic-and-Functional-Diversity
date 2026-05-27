generate_data <- function(p,
                          r,
                          signal = c("v1", "vp"),
                          family = c("gaussian", "poisson"),
                          c_val,
                          sigma2,
                          tree){
  
  signal <- match.arg(signal)
  family <- match.arg(family)
  
  Phy_DM <- getPhyloMatrix(tree, p)
  Phy_SM <- max(Phy_DM) - Phy_DM
  
  spec <- spectral_decomp(VC_phy_func = Phy_SM)
  V <- spec$P
  
  vec1 <- V[, 1]
  vecp <- V[, p]
  
  Lambda <- cbind(V[, 1], V[, 2])
  
  mu1 <- rep(0, p)
  
  mu2 <- if(signal == "v1"){
    c_val * vec1
  } else {
    c_val * vecp
  }
  
  Z1 <- matrix(rnorm(r * 2), r, 2)
  Z2 <- matrix(rnorm(r * 2), r, 2)
  
  U1 <- Z1 %*% t(Lambda)
  U2 <- Z2 %*% t(Lambda)
  
  Y1 <- matrix(0, r, p)
  Y2 <- matrix(0, r, p)
  
  for(j in 1:p){
    
    eta1 <- mu1[j] + U1[, j]
    eta2 <- mu2[j] + U2[, j]
    
    if(family == "gaussian"){
      Y1[, j] <- eta1 + rnorm(r, 0, sqrt(sigma2))
      Y2[, j] <- eta2 + rnorm(r, 0, sqrt(sigma2))
    }
    
    if(family == "poisson"){
      Y1[, j] <- rpois(r, lambda = exp(eta1))
      Y2[, j] <- rpois(r, lambda = exp(eta2))
    }
  }
  
  Y <- rbind(Y1, Y2)
  
  dat_long <- data.frame(
    id = rep(1:(2 * r), each = p),
    com = rep(rep(c(1, 2), each = r), each = p),
    sp = factor(rep(1:p, times = 2 * r)),
    y = c(t(Y))
  )
  
  dat_wide <- data.frame(
    com = factor(rep(c(1, 2), each = r)),
    Y
  )
  
  colnames(dat_wide)[-1] <- paste0("sp", 1:p)
  
  list(long = dat_long, wide = dat_wide)
}

fit_models <- function(dat, family = c("gaussian", "poisson"), globalTest, tree){
  
  family <- match.arg(family)
  
  dat_long <- dat$long
  
  fam <- if(family == "gaussian"){
    gaussian()
  } else if (family == "poisson") {
    nbinom2
  }
  
  if (globalTest == TRUE){
    fit_glmm_no_rr <- try(
      glmmTMB(
        y ~ com * sp,
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
    
    fit_glmm_no_rr_null <- try(
      glmmTMB(
        y ~ sp + com,
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
    
    fit_glmm_rr <- try(
      glmmTMB(
        y ~ com * sp +
          rr(sp + 0 | id, d = 2),
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
    
    fit_glmm_rr_null <- try(
      glmmTMB(
        y ~ sp + com +
          rr(sp + 0 | id, d = 2),
        data = dat_long,
        family = fam
      ),
      silent = TRUE
    )
  } else {
    p <- nlevels(dat$long$sp)
    Phy_DM <- getPhyloMatrix(tree, p)
    Phy_SM <- max(Phy_DM) - Phy_DM
    I_p <- diag(p) 
    J <- I_p - 1/p * matrix(1, p, p)   # centering matrix
    S_J <- t(J) %*% Phy_SM %*% J  # centered similarity matrix
    V_J <- spectral_decomp(VC_phy_func = S_J)$P  # eigenvectors of transformed matrix
    D <- spectral_decomp(VC_phy_func = S_J)$D
    val_num.eig <- get_num.eig(Methods = 'Method3', Methods_para = NA, D = D)
    fit_glmm_no_rr <- try(
      fit_model(type = 'p_mid', matrix_type = 'P_J', rr = FALSE,
                                 yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
                                 Distribution = fam, null = FALSE),
      silent = TRUE
    )
    
    fit_glmm_no_rr_null <- try(
      fit_model(type = 'p_mid', matrix_type = 'P_J', rr = FALSE,
                yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
                Distribution = fam, null = TRUE, com_left = NULL),
      silent = TRUE
    )
    
    fit_glmm_rr <- try(
      fit_model(type = 'p_mid', matrix_type = 'P_J', rr = TRUE,
                yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
                Distribution = fam, null = FALSE),
      silent = TRUE
    )
    
    fit_glmm_rr_null <- try(
      fit_model(type = 'p_mid', matrix_type = 'P_J', rr = TRUE,
                yX = dat_long, P_J = V_J, val_num.eig = val_num.eig, 
                Distribution = fam, null = TRUE, com_left = NULL),
      silent = TRUE
    )
  }
  
  get_lrt_p <- function(fit1, fit0){
    
    if(inherits(fit1, "try-error") ||
       inherits(fit0, "try-error")) return(NA)
    
    out <- try(anova(fit0, fit1), silent = TRUE)
    
    if(inherits(out, "try-error")) return(NA)
    
    out$`Pr(>Chisq)`[2]
  }
  
  data.frame(
    glmm_no_rr = get_lrt_p(fit_glmm_no_rr, fit_glmm_no_rr_null),
    glmm_rr    = get_lrt_p(fit_glmm_rr, fit_glmm_rr_null)
  )
}
