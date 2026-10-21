rm(list = ls())
load(file.path('data', "Para_PCV.RData"))

library(ape)
library(MASS)
source('R/MeanAbund.R')
source('R/U_current.R')
p <- 5
subsample <- sample(TruePara, p)
subtree <- keep.tip(tree, names(subsample))
DM <- cophenetic.phylo(subtree)
quantile = 0.35

r <- 64
alpha <- log(9/4)
beta <- log(4)
Eigen <- 3
NOS <- 1
FactorLoading <- Lambda

MeanAbund_Result <- MeanAbund(Subsample = subsample, DM, quantile, p,
  r, alpha, beta, Eigen, 
  FactorLoading, NOS)
Rao <- MeanAbund_Result$Rao
print(Rao - Rao[1])
MeanMatrix <- MeanAbund_Result$MeanMatrix
sp_Target <- MeanAbund_Result$sp_Target
sp<- colnames(DM)
length(unique(sp[sp_Target]))
plot(subtree)
