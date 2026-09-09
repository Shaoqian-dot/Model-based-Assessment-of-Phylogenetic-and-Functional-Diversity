rm(list = ls())
library(mvtnorm)
source('R/simulate_random_effects.R')
source('R/Boxplot.R')
source('R/Biplot.R')
set.seed(123)
logit <- function (x){
  fun <- log(x / (1 - x))
  return(fun)
}
expit <- function (x){
  fun <- exp(x)/ (1 + exp(x))
  return(fun)
}

# Design
r <- 128 # sample size in each com
p <- 5 # the number of species
q <- 6 # the number of coms
## Mean abundance
mu=matrix(c(0.8,0.2,0.8,0.2,0.8,0.2,0.8,0.8,0.2,0.8,0.8,0.8,0.2,0.2,0.8,0.8,0.8,0.8,0.2,0.2,2/3,1/9,2/3,1/9,2/3, 0.8, 0.2, 0.2, 0.8, 0.8),ncol=6)
rownames(mu)=LETTERS[1:5]
colnames(mu)=1:q
mu <- t(mu)
mu_site <- mu[rep(c(1 : q), each = r), ]

# Correlation
D=matrix(c(0,1,2,3,4,1,0,2,3,4,2,2,0,3,4,3,3,3,0,4,4,4,4,4,0),ncol=5)
colnames(D)=rownames(D)=LETTERS[1:5]
S <- max(D) - D 
mean <- rep(0, p)
sigma <- S
# 生成 r * p 个多元正态随机样本
resi_corr <- rmvnorm(n = r * q, mean = mean, sigma = sigma)
eta <- logit(mu_site) + resi_corr 
mu_eta <- expit(eta)

## Data generation
abundance <- mu_eta
abundance[] <- rbinom(length(mu_eta), size = 1, prob = mu_eta)

# Create community factor (1:q repeated r times each)
com <- factor(rep(1:q, each=r))

# Combine abundance matrix and community labels into a data frame
yX <- data.frame(abundance, com)

# Reshape data frame from wide to long format:
#   y: abundance value
#   sp: species
yX <- reshape(yX, direction = "long", varying = colnames(abundance),
              v.names = "y", timevar = "sp")

# Convert species column to factor
yX$sp <- factor(yX$sp)
 
## Model 
model <- glmmTMB (y ~ com * sp + rr(sp + 0 | id, 2),
                  family = 'binomial',
                  data = yX) 
## Diversity 
eta_T_fe_hat <- predict(model, 
                   type = "link",
                   re.form = NA)
eta_T_fe_hat_M <- matrix(
  as.numeric(eta_T_fe_hat),
  ncol = p,
  byrow = FALSE
)
eta_R_fe_hat <-  eta_T_fe_hat - apply(eta_T_fe_hat_M, 1, mean) 

RR <- simulate_random_effects(model = model)

eta_R_hat <- eta_R_fe_hat + RR

mu_eta_R_hat <- expit(eta_R_hat)

Div_Model <-  rowSums(mu_eta_R_hat %*% D * mu_eta_R_hat)
Div_Model_data <- data.frame(
  group = com,
  value = Div_Model
)
## True Abundance
rowSums(mu %*% D * mu)
## Boxplot
Boxplot(data = Div_Model_data,
        group_levels = c(
          '1',
          '2',
          '3',
          '4',
          '5',
          '6'),
        ylab = "mu_B^{T}Dmu_B"
        )
## Biplot 

Biplot(abundance = cbind(com, eta_R_hat), P = eigen(S)$vectors)
Biplot(abundance = cbind(com, expit(eta_R_hat)), P = eigen(S)$vectors)

J <- diag(rep(1, p)) -  1/p * matrix(1, p, p)
Biplot(abundance = cbind(com, eta_R_hat), P = eigen(S)$vectors)
Biplot(abundance = cbind(com, eta_R_hat), P = eigen(J %*% S %*% J)$vectors)
Biplot(abundance = cbind(com, expit(eta_R_hat)), P = eigen(J %*% S %*% J)$vectors)

