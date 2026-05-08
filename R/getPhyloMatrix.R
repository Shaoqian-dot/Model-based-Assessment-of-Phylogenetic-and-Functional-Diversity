# Function to construct a phylogenetic distance matrix for a random subset of species
# tree: a phylogenetic tree object (of class 'phylo')
# m: number of species to randomly select
getPhyloMatrix <- function(tree, m) {
  
  # Extract all species (tip labels) from the phylogenetic tree
  species <- tree$tip.label
  
  # Randomly select m species without replacement
  selected_species <- sample(species, m, replace = FALSE)
  
  # Prune the tree to keep only the selected species
  subtree <- keep.tip(tree, selected_species)
  
  # Compute the pairwise phylogenetic distance matrix (cophenetic distances)
  DM_phy_func <- cophenetic.phylo(subtree)
  
  # Return both the distance matrix
  return(DM_phy_func)
}