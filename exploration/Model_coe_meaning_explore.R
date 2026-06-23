library(glmmTMB)

set.seed(123)

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

# -----------------------------
# 2. 拟合 glmmTMB 模型
# -------------------- ---------

fit <- glmmTMB(
  y ~ sp + com : sp,
  data = dat,
  family = gaussian()
)
summary(fit)

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
