# Function to construct a phylogenetic distance matrix for a random subset of species
# tree: a phylogenetic tree object (of class 'phylo')
# m: number of species to randomly select
getPhyloMatrix <- function(tree, m) {
  
  # Extract all species (tip labels) from the phylogenetic tree
  species <- tree$tip.label
  
  # Randomly select m species without replacement
  selected_species <- sample(species, m, replace = FALSE)
  
  # Prune the tree to keep only the selected species
  subtree <<- keep.tip(tree, selected_species)
  
  # Compute the pairwise phylogenetic distance matrix (cophenetic distances)
  DM_phy_func <- cophenetic.phylo(subtree)
  
  # Return both the distance matrix
  return(DM_phy_func)
}

# Function to determine the number of eigenvalues to retain based on different methods
# Methods: a vector of method names ('Method1', 'Method2', ..., 'Method5')
# Methods_para: a vector of parameters corresponding to each method
# D: a square matrix (assumed to be diagonal with eigenvalues on the diagonal)
get_num.eig <- function (Methods, Methods_para, D){
  
  # Initialize a vector to store the number of eigenvalues for each method
  num.eig <- c()
  
  # Number of eigenvalues (assumed to be number of rows of D)
  m <- nrow(D)
  
  # Extract eigenvalues from the diagonal of D
  eig_vals <- diag(D)
  
  # Ensure Methods and Methods_para have the same length
  if (length(Methods_para) != length(Methods)) {
    stop("Methods_para must have the same length as Methods")
  }
  
  # Loop over each method
  for (i in 1 : length(Methods)){
    
    # Calculate proportion of variance explained (PVE) for each eigenvalue
    # PVE[j] = sum of first j eigenvalues / total sum of eigenvalues
    PVE <- numeric(m) 
    for (j in 1:m){
      PVE[j] <- sum(eig_vals[1:j]) / sum(eig_vals)
    }
    
    # Index vector and ratio index for cumulative comparison
    Index <- c(1 : m)
    Ratio_index <- Index/m
    
    # Apply the method to determine the number of eigenvalues
    if (Methods[i] == 'Method1'){
      # Method1: Use the provided parameter directly as the number of eigenvalues
      num.eig[i] <- Methods_para[i]
      
    } else if (Methods[i] == 'Method2'){
      # Method2: Retain enough eigenvalues to reach a specified proportion of variance
      val_num.eig <- 1
      while((sum(eig_vals[1 : val_num.eig]) / sum(eig_vals)) < Methods_para[i]){
        val_num.eig <- val_num.eig + 1
      }
      num.eig[i] <- val_num.eig
      
    } else if (Methods[i] == 'Method3'){
      # Method3: Retain eigenvalues greater than the mean eigenvalue
      num.eig[i] <- length(which(eig_vals > mean(eig_vals)))
      
    } else if (Methods[i] == 'Method4'){
      # Method4: Retain eigenvalues until max(PVE, cumulative index/m) exceeds 0.5
      num.eig[i] <- which(pmax(PVE, Ratio_index) > 1/2)[1]
      
    } else if (Methods[i] == 'Method5'){
      # Method5: Retain eigenvalues until min(PVE, cumulative index/m) exceeds 0.5
      num.eig[i] <- which(pmin(PVE, Ratio_index) > 1/2)[1]
    }
  }
  
  # Return the vector of number of eigenvalues for each method
  return(num.eig)
}

library(ape)
#source("R/getPhyloMatrix.R")
p <- 5
D <- matrix(c(0, 1, 2, 4, 4,
              1, 0, 2, 4, 4,
              2, 2, 0, 4, 4,
              4, 4, 4, 0, 2,
              4, 4, 4, 2, 0), 5, 5)

# An real example distance matrix
############################################################
### Load phylogenetic tree
############################################################
data_dir <- "data"

tree_file <- file.path(data_dir, "example.tre")

if(!file.exists(tree_file)){
  stop("Tree file not found: ", tree_file)
}

tree <- ape::read.tree(tree_file)

D <- getPhyloMatrix(tree = tree, m = p)
plot(subtree)

#S <- 1 - D/(max(D) + 1)
S <- max(D) - D
eigen <- eigen(S)
D_eigen <- diag(eigen$values)
D_eigen
V <- eigen$vectors
J <- diag(rep(1, p)) - 1/p * matrix(1, p, p)
V_sub <- J %*% V
qr(V_sub)$rank

S_new <- J %*% S %*% J
eigen_new <- eigen(S_new)
V_new <- eigen_new$vectors

# 求矩阵V_sub的一个极大线性无关组
q <- qr(V_sub)
pivot_cols <- q$pivot[1:q$rank]
V
V_new
V_sub
V_sub[, pivot_cols]

# Calculate the number of eligible contrasts under Method 3
num <- get_num.eig("Method3", NA, D_eigen)
num
D_eigen_new <- diag(eigen_new$values)
D_eigen_new
num_new <- get_num.eig("Method3", NA, D_eigen_new)
num_new
