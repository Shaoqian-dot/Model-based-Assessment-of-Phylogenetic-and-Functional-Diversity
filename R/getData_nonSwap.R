generate_data <- function(p,
                          r,
                          signal = c("v1", "vp"),
                          family = c("gaussian", "poisson"),
                          c_val,
                          sigma2,
                          tree,
                          V){ # Eigenvector matrix of the similarity matrix
  
  signal <- match.arg(signal)
  family <- match.arg(family)
  
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