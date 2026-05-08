# Function to produce u_ij
U <- function(Model, r, q, m){
  U <- matrix(NA, r * q, m)
  {
    if (sum(is.na(Model)) == 1 & length(is.na(Model)) == 1){
      mu_hat <- c(-0.21, 0.42)
      Sigma_hat <- matrix(c(0.61, 0.06, 
                            0.06, 0.99), 2, 2)
    }
    else{
      rep <- Model$obj$env$report()
      Lambda <- t(rep$fact_load[[3]])
      mu_hat <- rowMeans(Lambda)
      Sigma_hat <- cov(t(Lambda))
    }
  }
    Lambda_sample <- mvrnorm(m, mu_hat, Sigma_hat)
  for (i in 1 : (r * q)){
    for ( j in 1 : m){
      U[i, j] <- mvrnorm(1, rep(0, 2), diag(c(1, 1))) %*%
        Lambda_sample[j, ]
    }
  }
  return(U)
}
