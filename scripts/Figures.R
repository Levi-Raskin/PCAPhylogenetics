library(ggplot2)
library(ggtree)

output <- "/Users/levir/Documents/GitHub/PCAPhylogenetics/manuscript/figures"

treeSubset <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/sampledTrees.tsv")

#convert to phylo object
phyloList <- list()
for(i in 1:nrow(treeSubset)){
  phyloList[[i]] <- ape::read.tree(text = as.character(treeSubset[i, 2]))
}


# Figure: BM continuous traits --------------------------------------------

tree <- phyloList[[1]]
tree$tip.label <- gsub("_", tree$tip.label, replacement = " ")
tree <- ggtree()