Biplot <- function(abundance, P) {
     
       # ---------------------------------------------------------------------------
     # Extract species-level effects
       # ---------------------------------------------------------------------------
     
       eta_M <- as.matrix(
           abundance[, setdiff(colnames(abundance), "com")]
         )
       
         # Ensure that the matrix is numeric.
         eta_M <- matrix(
             as.numeric(eta_M),
             nrow = nrow(eta_M),
             ncol = ncol(eta_M),
             dimnames = dimnames(eta_M)
           )
         
           
           # ---------------------------------------------------------------------------
         # Project species-level effects onto the first two eigenvectors
           # ---------------------------------------------------------------------------
         aim <- eta_M %*% P[, c(1, 2)]
         
           colnames(aim) <- paste0(
               "v",
               seq_len(ncol(aim))
             )
           
             
             # ---------------------------------------------------------------------------
           # Create plotting data
             # ---------------------------------------------------------------------------
           
             df_scores <- data.frame(
                 com = abundance[, "com"],
                 score1 = aim[, 1],
                 score2 = aim[, 2]
               )
             
               df_scores$com <- factor(
                   df_scores$com
                 )
               
                 
                 # ---------------------------------------------------------------------------
               # Plot
                 # ---------------------------------------------------------------------------
               
                 ggplot(
                     df_scores,
                     aes(
                         x = score1,
                         y = score2,
                         color = com
                       )
                   ) +
                   geom_point(size = 3) +
                   labs(
                       x = expression(v[1]^T * eta),
                       y = expression(v[2]^T * eta)
                     ) +
                   coord_fixed(ratio = 1) +
                   theme_classic()
             }