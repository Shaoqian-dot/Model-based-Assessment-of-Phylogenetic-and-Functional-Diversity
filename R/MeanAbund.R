# The function swaps multiple pairs of mean abundances; 
# The actual number of swaps depends on how many indexes in Index_No_change 
# and how many pairs of species (spRef_Small and sp_Large) there are outside of Index_No_change.

MeanAbund <- function (Subsample, DM, quantile, p,
                       r, alpha, beta, Eigen, 
                       FactorLoading, NOS, V_J, family_type) { # NOS : number of swaps
 
#################################################### Experiment Parameters
  # TruePara <- mean_abundance
  # quantile <- 0.35
  # p <- 5
  # alpha <- 2.426285
  # beta <- 2.426285
  # model <- model_updated
####################################################  
 
  #family_type <- "binomial"  # 或 "poisson"
  
  family_obj <- switch(
    family_type,
    binomial = stats::binomial(link = "logit"),
    poisson  = stats::poisson(link = "log"),
    stop("Unsupported family_type: ", family_type)
  )
  
  link_fun <- family_obj$linkfun
  inv_link <- family_obj$linkinv
  
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
 NOS <- min(length(Index_No_change), NOS) # The actual NOS
 sp_Target <- c()
 for (i in 1 : NOS) {
   
   ## Another way to choose species with small and large distances.
   ## Randomly choose a pair of phylogenetically/functionally equal species
   Target <- sample(1 : length(Index_No_change), 1) 
   Index_Target <- Index_No_change[[Target]] 
   spRef <- sp[Index_Target][1]
   spRef_equal <- sp[Index_Target][2]
   
   ## Randomly choose two species which are phylogenetically/functionally close and far from spRef
   DisSort <- sort(DM[spRef,
                      -c(unlist(Index_No_change, use.names = FALSE), sp_Target)]) # Choose spRef_Small and spRef_Large outside of Index_No_change
   
   spOrder <- names(DisSort)
   spRef_Small <- spOrder[sample(1 : round(length(DisSort) * quantile), 1)]
   spRef_Large <- spOrder[sample((round(length(DisSort) * quantile) + 1) : length(DisSort), 1)]
   
   if (length(unique(c(spRef, spRef_equal, spRef_Small, spRef_Large))) < 4 ){
     DisSort <- sort(DM[spRef,-c(Index_Target, sp_Target)])
     spOrder <- names(DisSort)
     spRef_Small <- spOrder[sample(1 : round(length(DisSort) * quantile), 1)]
     spRef_Large <- spOrder[sample((round(length(DisSort) * quantile) + 1) : length(DisSort), 1)]
   }
     
   sp_Target <- c(sp_Target, match(c(spRef, spRef_equal, spRef_Small, spRef_Large), sp))
   if (i == 1) {
     # Construct Mean Matrix
     mu1 <- Subsample[sp]
     mu2 <- mu1
     mu3 <- mu1
     mu4 <- mu1
   } else {
     
   }
   
   Effect_low  <- inv_link(link_fun(0.1848) - (beta/2)) # 0.1848 is the overall mean abundance across the 115 species
   Effect_high <- inv_link(link_fun(0.1848) + (beta/2))
   
   mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_low,  Effect_high,
                                                             Effect_high,  Effect_high)
   
   mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high, Effect_low,
                                                             Effect_high,  Effect_high)
   
   mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
                                                             Effect_low,  Effect_high)
   
   mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
                                                             Effect_high, Effect_low)
   
   ref_indices <- match(c(spRef, spRef_equal, spRef_Small, spRef_Large), sp)
   
   # Note: Sometimes, the number of pairs chosen in Index_No_change is unequal to NOS because there are more than 2 phylogenetically equal species.  
   Index_No_change <- purrr::discard(
     Index_No_change,
     ~ any(.x %in% ref_indices)
   ) 
   
   # # eta_Mean <- mean(qlogis(mu1))
   # # 
   # # mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean - (beta/2)), plogis(eta_Mean + (beta/2)),
   # #                                                  plogis(eta_Mean + (beta/2)), plogis(eta_Mean + (beta/2)))
   # # 
   # # mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean + beta/2), plogis(eta_Mean - beta/2),
   # #                                                  plogis(eta_Mean + beta/2), plogis(eta_Mean + beta/2))
   # # 
   # # mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean + beta/2), plogis(eta_Mean + beta/2),
   # #                                                  plogis(eta_Mean - beta/2), plogis(eta_Mean + beta/2))
   # # 
   # # mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(plogis(eta_Mean + beta/2), plogis(eta_Mean + beta/2),
   # #                                                  plogis(eta_Mean + beta/2), plogis(eta_Mean - beta/2))
   # 
   # # mu_Mean <- mean(mu1)
   # # Effect_low <- plogis(qlogis(mu_Mean) - (beta/2))
   # # Effect_high <- plogis(qlogis(mu_Mean) + (beta/2))
   # # 
   # # mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_low,  Effect_high,
   # #                                                           Effect_high,  Effect_high)
   # # 
   # # mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high, Effect_low,
   # #                                                            Effect_high,  Effect_high)
   # # 
   # # mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high,  Effect_high,
   # #                                                           Effect_low,  Effect_high)
   # # 
   # # mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
   # #                                                            Effect_high, Effect_low)
   # 
   # # Effect_low <- mu1[spRef]
   # # Effect_high <- plogis(qlogis(mu1[spRef]) + beta)
   # # 
   # # mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_low,  Effect_high,
   # #                                                           Effect_high,  Effect_high)
   # # 
   # # mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high, Effect_low,
   # #                                                            Effect_high,  Effect_high)
   # # 
   # # mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c( Effect_high,  Effect_high,
   # #                                                            Effect_low,  Effect_high)
   # # 
   # # mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
   # #                                                           Effect_high, Effect_low)
   # 
   # Effect_low  <- plogis(qlogis(0.1848) - (beta/2)) # 0.1848 is the overall mean abundance across the 115 species
   # Effect_high <- plogis(qlogis(0.1848) + (beta/2))
   # 
   # mu1[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_low,  Effect_high,
   #                                                           Effect_high,  Effect_high)
   # 
   # mu2[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high, Effect_low,
   #                                                           Effect_high,  Effect_high)
   # 
   # mu3[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
   #                                                           Effect_low,  Effect_high)
   # 
   # mu4[c(spRef, spRef_equal, spRef_Small, spRef_Large)] <- c(Effect_high,  Effect_high,
   #                                                           Effect_high, Effect_low)
   # 
   # ref_indices <- match(c(spRef, spRef_equal, spRef_Small, spRef_Large), sp)
   # 
   # Index_No_change <- purrr::discard(
   #   Index_No_change,
   #   ~ any(.x %in% ref_indices)
   # )
 }
 
 mu5 <- inv_link(link_fun(mu1) + alpha)
 
 MeanMatrix <- rbind(mu1, mu2, mu3, mu4, mu5)
 rownames(MeanMatrix) <- paste0('com', 1 : nrow(MeanMatrix))
 
 # Rao <- rowSums(MeanMatrix %*% DM * MeanMatrix)
 # Correlation
 S <- max(DM) - DM
 J <- diag(1, p) - 1/p * matrix(1, p, p)
 U_corr <- U(p = p, r = r, q = nrow(MeanMatrix), Eigen = Eigen,
             FactorLoading = FactorLoading, sp_Sim = sp, V_J = V_J)
 
 # Full Mean Matrix
 stopifnot(identical(colnames(MeanMatrix), colnames(U_corr)))
 
 MeanMatrix_full <- inv_link(link_fun(MeanMatrix[rep(1: nrow(MeanMatrix), each = r), ]) + U_corr)
 #return(MeanMatrix_full)
 return(list(MeanMatrix_full = MeanMatrix_full,
             MeanMatrix = MeanMatrix))

 }


