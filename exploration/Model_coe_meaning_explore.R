library(glmmTMB)

#set.seed(123)

# -----------------------------
# 1. 生成模拟数据
# -----------------------------

n <- 300

dat <- data.frame(
  sp  = factor(sample(c("A", "B", "C"), n, replace = TRUE)),
  com = factor(sample(c("X", "Y"), n, replace = TRUE))
)

# 设置真实参数
beta0 <- 10

# sp effect
beta_sp <- c(
  A = 0,
  B = 2,
  C = -1
)

# com effect
beta_com <- c(
  X = 0,
  Y = 3
)

# interaction effect
beta_int <- matrix(
  c(
    0,  0,   # A
    0,  2,   # B
    0, -2    # C
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("A", "B", "C"),
                  c("X", "Y"))
)

# 构造均值
mu <- numeric(n)

for(i in 1:n){
  
  sp_i  <- dat$sp[i]
  com_i <- dat$com[i]
  
  mu[i] <-
    beta0 +
    beta_sp[sp_i] +
    beta_com[com_i] +
    beta_int[sp_i, com_i]
}

# 生成 normal response
sigma <- 2

dat$y <- rnorm(n, mean = mu, sd = sigma)
dat$tmp <- rnorm(n, mean = mu, sd = sigma)
# 查看前几行
head(dat)

DM <- matrix(c(0, 1, 2, 
               1, 0, 2,
               2, 2, 0), 3, 3)
SM <- 1 - DM/2

J <- diag(c(1, 1, 1)) - 1/3 * matrix(1, 3, 3)
SM_J <- J %*% SM %*% J
V_J <- eigen(SM_J)$vectors
# V_J[, c(1 : (p-1))] as contrasts, satisfying sum-to-zero
contrasts(dat$sp) = V_J[, c(1, 2)]

# Sum-to-zero contrast
contrasts(dat$sp) = contr.sum(3)

p <- 3
val_num.eig <- 1
Sp <- model.matrix(~0 + sp, data = dat)
Lsp <- Sp %*% V_J
colnames(Lsp) <- paste0("p", 1:3)

yX_tmp <- cbind(dat, Lsp)

K <- diag(p - val_num.eig + 1)
colnames(K) <- c("(Intercept)", paste0("p", (val_num.eig + 1) : p))
# rownames(K) <- paste0("p", (val_num.eig + 1) : p)

# -----------------------------
# 2. 拟合 glmmTMB 模型
# -------------------- ---------
# fit <- glmmTMB(
#   y ~ sp + com : tmp,
#   data = dat,
#   family = gaussian()
# )
# summary(fit)

fit_redu <- glmmTMB(
  y ~ p1 + com : p1 + diag(p2 | com),
  data = yX_tmp,
  family = gaussian()
)
summary(fit_redu)


fit_redu_com <- glmmTMB(
  y ~ p1 + com : p1 + diag(p2 | com),
  data = yX_tmp,
  family = gaussian()
)
summary(fit_redu_com)


fit_full <- glmmTMB(
  y ~ p1 + com : p1 + propto(p2 + p3| com, K),
  data = yX_tmp,
  family = gaussian()
)

summary(fit_full)


fit_full_com <- glmmTMB(
  y ~ com * p1 + propto(0 + p2 + p3 | com, K),
  data = yX_tmp,
  family = gaussian()
)

summary(fit_full_com)

fit <- glmmTMB(
  y ~ sp * com,
  data = dat,
  family = gaussian()
)
summary(fit)

fit_no_Intercept <- glmmTMB(
  y ~ 0 + sp * com,
  data = dat,
  family = gaussian()
)

summary(fit_no_Intercept)

############################################################ Use my simulation data set
library(ape)
library(MASS)
library(glmmTMB)
source("R/getPhyloMatrix.R")
source("R/spectral_decomp.R")
source("R/getData_current.R")
source("R/GroupingPairs.R")
source("R/mean_matrix_current.R")
source("R/sample_pairs.R")
source("R/U_current.R")
source("R/fit_model_current2.R")

p <- 5
NOPS <- 1
q <- 5
alpha <- log(9/4)
beta <- log(2)
r <- 64
quantile <- 0.35
Corr <- 1
Eigen <- 1
family_type <- 'poisson'
family_fit <- 'nbinom2'
val_num.eig <- 1
############################################################
### Load phylogenetic tree
############################################################

data_dir <- "data"

tree_file <- file.path(data_dir, "example.tre")

if(!file.exists(tree_file)){
  stop("Tree file not found: ", tree_file)
}

tree <- ape::read.tree(tree_file)

DM_phy_func <- getPhyloMatrix(tree = tree, m = p)
VC_phy_func <- 1 - DM_phy_func / max(DM_phy_func)

I_p <- diag(p)
J <- I_p - 1 / p * matrix(1, p, p)
S_J <- t(J) %*% VC_phy_func %*% J

eig_J <- spectral_decomp(VC_phy_func = S_J)
V_J <- eig_J$P

eig <- spectral_decomp(VC_phy_func = VC_phy_func)
V <- eig$P

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

Sp <- model.matrix(~0 + sp, data = dat_long)
Lsp <- Sp %*% V_J
colnames(Lsp) <- paste0("p", 1:p)

yX_tmp <- cbind(dat_long, Lsp)
yX_tmp$p6 <- 1
yX_tmp$p7 <- yX_tmp$p2 + yX_tmp$p3

K1 <- diag(p + 1  - val_num.eig)
K2 <- diag(p  - val_num.eig)

colnames(K1) <- c("(Intercept)", paste0("p", (val_num.eig + 1) : p))
colnames(K2) <- paste0("p", (val_num.eig + 1) : p)


fit_glmm_no_rr_com <- glmmTMB(y ~  com * p1 + propto(0 + p2 + p3 + p4 + p5 | com, K2),
                              data = yX_tmp,
                              family = family_fit
)
summary(fit_glmm_no_rr_com)

fit_glmm_no_rr1 <- glmmTMB(y ~ p1 + com:(p1) + propto(0 + p2 + p3 + p4 + p5 | com, K2),
                           data = yX_tmp,
                           family = family_fit
)
summary(fit_glmm_no_rr1)

fit_glmm_no_rr2 <- glmmTMB(y ~ (p1) + com:(p1) + propto(p2 + p3 + p4 + p5| com, K1),
                          data = yX_tmp,
                          family = family_fit
                          )
summary(fit_glmm_no_rr2)

fit_glmm_no_rr22 <- glmmTMB(y ~ (p1) + com:(p1) + propto(0 + p2 + p3 + p4 + p5 | com, K2),
                           data = yX_tmp,
                           family = family_fit
)
summary(fit_glmm_no_rr22)


fit_glmm_rr <- glmmTMB(y ~ p1 + com:(p1) + propto(0 + p2 + p3 + p4 + p5 | com, K2) +
                         rr(sp + 0|id, 2),
                           data = yX_tmp,
                           family = family_fit
)
summary(fit_glmm_rr)
