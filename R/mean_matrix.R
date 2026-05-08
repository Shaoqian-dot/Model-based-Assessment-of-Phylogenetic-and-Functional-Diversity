# Create a function for producing a mean matrix
mean_matrix <- function (Index_change, beta, alpha, 
                         q, m, r, NOPS, Model, Corr){
  {
    if (m %% 2 == 0){
      Means_logit <- c(rep(c(-beta, beta), times = (q-1)*(m/2)),
                       rep(c(alpha - beta, alpha + beta), times = m/2))
      
      Mean_matrix_logit <- matrix(Means_logit, c(q, m), byrow = TRUE)
      
    } else {
      
      Means_logit <- c(rep(c(-beta, beta), times = (q-1)*((m-1)/2)),
                       rep(c(alpha - beta, alpha + beta), times = (m-1)/2))
      
      Mean_matrix_logit <- matrix(Means_logit, c(q, (m-1)), byrow = TRUE)
      
      Mean_matrix_logit <- cbind(
        Mean_matrix_logit,
        c(rep(-beta, times = q-1), alpha - beta)
      )
    }
  }
  NOPS_rep <- NOPS # An intermediate variable
  for (i in 2 : (q - 1)){
    {# Reduce the number of swaps if there are not enough pairs for swapping (less than NOPS) in No PD group
      if ((length(Index_change[[i-1]])/2 < NOPS) & (length(Index_change[[i-1]]) != 0)){
        NOPS <- length(Index_change[[i-1]])/2
      }
      else if (length(Index_change[[i-1]])/2 >= NOPS) {
        
      }
      else  {
        message(sprintf("No satisfied pair in group %d", i - 1))
      }
    }
    #An alternative way of swapping means, but it doesn't work when there is NULL in Index_change
    #Index_aim <- Index_change[[i-1]][c(1:(NOPS*2))]
    #Mean_matri_temp <- Mean_matri[i, ]# A temporary vector
    #Mean_matri[i, Index_aim[seq(1, length(Index_aim), 2)]] <-
    #Mean_matri_temp[Index_aim[seq(2, length(Index_aim), 2)]]# NOPS is the number of species pairs we want to swap; 2: pair; The sum of two adjacent mean abundances in com1 should be 1.
    #Mean_matri[i, Index_aim[seq(2, length(Index_aim), 2)]] <-
    #Mean_matri_temp[Index_aim[seq(1, length(Index_aim), 2)]]# NOPS is the number of species pairs we want to swap; 2: pair; The sum of two adjacent mean abundances in com1 should be 1.
    #An alternative way of swapping means, but it doesn't work when there is NULL in Index_change
    Mean_matrix_logit[i, Index_change[[i-1]][c(1:(NOPS*2))]] <-
      -beta + beta - Mean_matrix_logit[i, Index_change[[i-1]][c(1: (NOPS*2))]]# NOPS is the number of species pairs we want to swap; 2: pair; The sum of two adjacent mean abundances in com1 should be 1.
    NOPS <- NOPS_rep
  }
  assign("Mean_matrix_logit", Mean_matrix_logit, envir = .GlobalEnv)# Set yX a global variable
  {
    if (Corr == 1){
      Model<- Model # The model is a model from real applications with diag(0 + com | sp)
      U_M <- U(Model, r, q, m)
      # Linear predictor
      Mean_matri_full <- plogis(Mean_matrix_logit[rep(1 : q, each = r), ] + U_M)
    }
    else {
      Mean_matri_full <- plogis(Mean_matrix_logit[rep(1 : q, each = r), ])
    }
  }
  assign('Mean_matri_full', Mean_matri_full, envir = .GlobalEnv)
  return(Mean_matri_full)
}
# Create a function for producing a mean matrix