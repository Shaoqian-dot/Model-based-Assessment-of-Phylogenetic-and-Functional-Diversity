MeanAbund <- function (Subsample, DM, quantile, p,
                       r, alpha, beta, Eigen, FactorLoading) {
 
#################################################### Experiment Parameters
  # TruePara <- mean_abundance
  # quantile <- 0.35
  # p <- 5
  # alpha <- 2.426285
  # beta <- 2.426285
  # model <- model_updated
####################################################  
 sp <- colnames(DM)
 
 # Randomly Choose Species to Swap
 ## Choose All the Pairs of Indexes Pointing to phylogenetically/functionally species
 Index_No_change <- list()
 k <- 1
 for (i in 1 : (p - 1)){
   for (j in (i + 1) : p) {
     if (isTRUE(all.equal(
       unname(DM[i, -c(i, j)]),
       unname(DM[j, -c(j, i)])
     ))){
       Index_No_change[[k]] <- c(i, j)
       k <- k + 1
     }
   }
 }
 ## Randomly choose a pair of phylogenetically/functionally equal species
 Index_Target <- Index_No_change[[sample(1 : (k-1), 1)]] # k - 1 is the number of pairs containing phylogenetically equal species 
 spRef <- sp[Index_Target][1]
 spRef_equal <- sp[Index_Target][2]
 
 ## Randomly choose two species which are phylogenetically/functionally close and far from spRef
 DisSort <- sort(DM[spRef,-Index_Target]) 
 spOrder <- names(DisSort)
 spRef_Small <- spOrder[sample(1 : round((p - 2) * quantile), 1)] # p - 2 is the length of DisSort
 spRef_Large <- spOrder[sample((round((p - 2) * quantile) + 1) : (p - 2), 1)]
 
 # Construct Mean Matrix
 mu1 <- Subsample[sp]
 mu2 <- mu1
 mu3 <- mu1
 mu4 <- mu1
 # eta_Mean <- mean(qlogis(mu1))
 # 
 # mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean - (beta/2)), plogis(eta_Mean + (beta/2)),
 #                                                  plogis(eta_Mean + (beta/2)), plogis(eta_Mean + (beta/2)))
 # 
 # mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean + beta/2), plogis(eta_Mean - beta/2),
 #                                                  plogis(eta_Mean + beta/2), plogis(eta_Mean + beta/2))
 # 
 # mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean + beta/2), plogis(eta_Mean + beta/2),
 #                                                  plogis(eta_Mean - beta/2), plogis(eta_Mean + beta/2))
 # 
 # mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean + beta/2), plogis(eta_Mean + beta/2),
 #                                                  plogis(eta_Mean + beta/2), plogis(eta_Mean - beta/2))
 
 # mu_Mean <- mean(mu1)
 # Effect_low <- plogis(qlogis(mu_Mean) - (beta/2))
 # Effect_high <- plogis(qlogis(mu_Mean) + (beta/2))
 # 
 # mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_low,  Effect_high,
 #                                                           Effect_high,  Effect_high)
 # 
 # mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high, Effect_low,
 #                                                            Effect_high,  Effect_high)
 # 
 # mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high,  Effect_high,
 #                                                           Effect_low,  Effect_high)
 # 
 # mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
 #                                                            Effect_high, Effect_low)
 
 Effect_low <- mu1[spRef]
 Effect_high <- plogis(qlogis(mu[spRef]) + beta)
 
 mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_low,  Effect_high,
                                                           Effect_high,  Effect_high)
 
 mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high, Effect_low,
                                                            Effect_high,  Effect_high)
 
 mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high,  Effect_high,
                                                            Effect_low,  Effect_high)
 
 mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
                                                           Effect_high, Effect_low)
 
 
 mu5 <- plogis(qlogis(mu1) + alpha)
 
 MeanMatrix <- rbind(mu1, mu2, mu3, mu4, mu5)
 rownames(MeanMatrix) <- paste0('com', 1 : nrow(MeanMatrix))
 
 # Correlation
 S <- max(DM) - DM
 J <- diag(1, p) - 1/p * matrix(1, p, p)
 U_corr <- U(p = p, r = r, q = nrow(MeanMatrix), Eigen = Eigen,
             FactorLoading = FactorLoading, sp_Sim = sp)
 
 # Full Mean Matrix
 stopifnot(identical(colnames(MeanMatrix), colnames(U_corr)))
 
 MeanMatrix_full <- plogis(qlogis(MeanMatrix[rep(1: nrow(MeanMatrix), each = r), ]) + U_corr)
 return(MeanMatrix_full)

 }


