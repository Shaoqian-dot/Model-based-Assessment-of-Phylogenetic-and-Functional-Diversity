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

S <- 1 - D/(max(D) + 1)
eigen <- eigen(S)
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

