# Setup -------------------------------------------------------------------
library(ape)
library(dplyr)
library(ggplot2)
library(ggtree)
library(phangorn)
library(phytools)
library(tictoc)
library(data.table)
library(geomorph)

setwd("/Users/levir/Documents/GitHub/PCAPhylogenetics")

set.seed(5)

rSPRFunc <- function(tree1, tree2){
  treelist <- list(tree1, tree2)
  write.tree(treelist, paste(getwd(), "/rSPR/res.nwk", sep = ""))
  res <- system(
    paste(getwd(), 
          "/rSPR/rspr -pairwise -unrooted -no-symmetric-pairwise < ",
          getwd(), 
          "/rSPR/res.nwk", sep = "") 
    , intern = T)
  nums <- sub("^[^,]*,\\s*([0-9]+)\"?$", "\\1", res[[1]])
  return(as.integer(nums))
}

rSPRFunc2 <- function(tree1, tree2){
  treelist <- list(tree1, tree2)
  write.tree(treelist, paste(getwd(), "/rSPR/res2.nwk", sep = ""))
  res <- system(
    paste(getwd(), 
          "/rSPR/rspr2 -pairwise -unrooted -no-symmetric-pairwise < ",
          getwd(), 
          "/rSPR/res2.nwk", sep = "") 
    , intern = T)
  nums <- sub("^[^,]*,\\s*([0-9]+)\"?$", "\\1", res[[1]])
  return(as.integer(nums))
}

rSPRFunc3 <- function(tree1, tree2){
  treelist <- list(tree1, tree2)
  write.tree(treelist, paste(getwd(), "/rSPR/res3.nwk", sep = ""))
  res <- system(
    paste(getwd(), 
          "/rSPR/rspr3 -pairwise -unrooted -no-symmetric-pairwise < ",
          getwd(), 
          "/rSPR/res3.nwk", sep = "") 
    , intern = T)
  nums <- sub("^[^,]*,\\s*([0-9]+)\"?$", "\\1", res[[1]])
  return(as.integer(nums))
}

rSPRFunc4 <- function(tree1, tree2){
  treelist <- list(tree1, tree2)
  write.tree(treelist, paste(getwd(), "/rSPR/res4.nwk", sep = ""))
  res <- system(
    paste(getwd(), 
          "/rSPR/rspr4 -pairwise -unrooted -no-symmetric-pairwise < ",
          getwd(), 
          "/rSPR/res4.nwk", sep = "") 
    , intern = T)
  nums <- sub("^[^,]*,\\s*([0-9]+)\"?$", "\\1", res[[1]])
  return(as.integer(nums))
}


sprMast <- function(tree1, tree2){
  sprDist <- rSPRFunc(tree1, tree2)
  mastSize <- Ntip(phangorn::mast(tree1, tree2))
  return(sprDist/mastSize)
}

sprMast2 <- function(tree1, tree2){
  sprDist <- rSPRFunc2(tree1, tree2)
  mastSize <- Ntip(phangorn::mast(tree1, tree2))
  return(sprDist/mastSize)
}


#read tree trace
# trees <- data.table::fread("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/hominin.trees",
                           # nThread = 10)
#burn in
# trees <- trees[round(0.1 * nrow(trees)) : nrow(trees),]
# sample 1,000 trees from the posterior distribution
# treeSubset <- trees[sample(1:nrow(trees), 1000, replace = FALSE), ]
# rm(trees)
# write.csv(treeSubset,"/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/sampledTrees.csv")
treeSubset <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/sampledTrees.tsv")

#convert to phylo object
phyloList <- list()
for(i in 1:nrow(treeSubset)){
  phyloList[[i]] <- ape::read.tree(text = as.character(treeSubset[i, 2]))
}

# Parameters --------------------------------------------------------------

numDatasets <- 1:25
numCharacters <- c(10, 100, 250)
setRate <- c(0.1, 1, 10)
proportionConflicting <- c(0, 0.1, 0.25, 0.5) #proportion of characters that exhibit conflicting phylogenetic signal
# proportionConflicting <- c(0)
numTaxaConflicting <- c(2) #just rearranging 2 taxa to start
#variableRates # drawn from gamma with parameters (1,10); (1,1); (10,1)

# # Simulation --------------------------------------------------------------
# 
# simulationFunction <- function(tree, treeIdx){
#   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
#   
#   for(nC in numCharacters){
#     for(sR in setRate){
#       for(pc in proportionConflicting){
#         for(d in numDatasets){
#           
#           datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#           rownames(datasetD) = tree$tip.label
#           for(t in 1:nC){
#             traits <- fastBM(tree, sig2 = sR)
#             datasetD[,t] <- traits
#           }
#           
#           whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#           for(t in whichConflicting){
#             newVec <- datasetD[,t]
#             tax <- sample(nrow(datasetD), 2, replace = F)
#             tax1Val <- newVec[tax[1]]
#             tax2Val <- newVec[tax[2]]
#             
#             datasetD[tax[1],t] <- tax2Val
#             datasetD[tax[2],t] <- tax1Val
#           }
#           dt <- as.data.table(datasetD)
#           rownames(dt) <- rownames(datasetD)
#           
#           data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#           pca <- prcomp(datasetD)
#           plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#           pcDistMat<- dist(plotdf)
#           njTree <- nj(pcDistMat)
#           unrootedTree <- unroot(tree)
#           unrootedNJTree <- unroot(njTree)
#           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
#           names(resVec) <- NULL
#           resMat[nrow(resMat) + 1, ] <- resVec
#         }
#       }
#     }
#     for(vR in 1:3){
#       for(pc in proportionConflicting){
#         if(vR == 1){
#           # 1, 10
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 1, scale = 10)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else if (vR == 2){
#           # 1, 1
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 1, scale = 1)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else{
#           # 10, 1
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 10, scale = 1)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             data.table::fwrite((dt),  paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }
#       }
#     }
#   }
#   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# }
# 
# resList <- list()
# 
# for(i in 1:length(phyloList)){
#   print(paste("On tree idx: ", i))
#   resList[[i]] <- simulationFunction(phyloList[[i]], i)
# }
# 
# saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults.rds")
# 
# simulationFunctionAllPCs <- function(tree, treeIdx){
#   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
#   
#   for(nC in numCharacters){
#     for(sR in setRate){
#       for(pc in proportionConflicting){
#         for(d in numDatasets){
#           datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#           pca <- prcomp(datasetD)
#           pcDistMat<- dist(pca$x)
#           njTree <- nj(pcDistMat)
#           unrootedTree <- unroot(tree)
#           unrootedNJTree <- unroot(njTree)
#           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
#           names(resVec) <- NULL
#           resMat[nrow(resMat) + 1, ] <- resVec
#         }
#       }
#     }
#     for(vR in 1:3){
#       for(pc in proportionConflicting){
#         if(vR == 1){
#           # 1, 10
#           for(d in numDatasets){
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else if (vR == 2){
#           # 1, 1
#           for(d in numDatasets){
#             
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else{
#           # 10, 1
#           for(d in numDatasets){
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }
#       }
#     }
#   }
#   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# }
# 
# resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
#   return(simulationFunctionAllPCs(phyloList[[i]], i))
# })
# 
# saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsAllPCs.rds")
# 
# # Simulation-- 25, 50 characters -----------------------------------------
# 
# numCharacters <- c(25, 50)
# 
# simulationFunction <- function(tree, treeIdx){
#   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
#   
#   for(nC in numCharacters){
#     for(sR in setRate){
#       for(pc in proportionConflicting){
#         for(d in numDatasets){
#           
#           datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#           rownames(datasetD) = tree$tip.label
#           for(t in 1:nC){
#             traits <- fastBM(tree, sig2 = sR)
#             datasetD[,t] <- traits
#           }
#           
#           whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#           for(t in whichConflicting){
#             newVec <- datasetD[,t]
#             tax <- sample(nrow(datasetD), 2, replace = F)
#             tax1Val <- newVec[tax[1]]
#             tax2Val <- newVec[tax[2]]
#             
#             datasetD[tax[1],t] <- tax2Val
#             datasetD[tax[2],t] <- tax1Val
#           }
#           
#           dt <- as.data.table(datasetD)
#           rownames(dt) <- rownames(datasetD)
#           data.table::fwrite((dt), paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#           pca <- prcomp(datasetD)
#           plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#           pcDistMat<- dist(plotdf)
#           njTree <- nj(pcDistMat)
#           unrootedTree <- unroot(tree)
#           unrootedNJTree <- unroot(njTree)
#           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
#           names(resVec) <- NULL
#           resMat[nrow(resMat) + 1, ] <- resVec
#         }
#       }
#     }
#     for(vR in 1:3){
#       for(pc in proportionConflicting){
#         if(vR == 1){
#           # 1, 10
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 1, scale = 10)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else if (vR == 2){
#           # 1, 1
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 1, scale = 1)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else{
#           # 10, 1
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 10, scale = 1)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt,  paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }
#       }
#     }
#   }
#   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# }
# 
# resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
#   return(simulationFunction(phyloList[[i]], i))
# })
# 
# saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults25char50char.rds")
# 
# simulationFunctionAllPCs <- function(tree, treeIdx){
#   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
#   
#   for(nC in numCharacters){
#     for(sR in setRate){
#       for(pc in proportionConflicting){
#         for(d in numDatasets){
#           datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#           pca <- prcomp(datasetD)
#           pcDistMat<- dist(pca$x)
#           njTree <- nj(pcDistMat)
#           unrootedTree <- unroot(tree)
#           unrootedNJTree <- unroot(njTree)
#           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
#           names(resVec) <- NULL
#           resMat[nrow(resMat) + 1, ] <- resVec
#         }
#       }
#     }
#     for(vR in 1:3){
#       for(pc in proportionConflicting){
#         if(vR == 1){
#           # 1, 10
#           for(d in numDatasets){
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else if (vR == 2){
#           # 1, 1
#           for(d in numDatasets){
#             
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else{
#           # 10, 1
#           for(d in numDatasets){
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }
#       }
#     }
#   }
#   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# }
# 
# resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
#   return(simulationFunctionAllPCs(phyloList[[i]], i))
# })
# 
# saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsAllPCs25char50char.rds")
# 
# # mahalanobis_dist_matrix <- function(data) {
# #   # Ensure data is a matrix
# #   data <- as.matrix(data)
# #   
# #   # Compute the covariance matrix of the data
# #   cov_matrix <- cov(data)
# #   
# #   # Initialize distance matrix
# #   n <- nrow(data)
# #   dist_matrix <- matrix(0, n, n)
# #   rownames(dist_matrix) <- rownames(data)
# #   colnames(dist_matrix) <- rownames(data)
# #   
# #   # Compute pairwise Mahalanobis distances
# #   for (i in 1:n) {
# #     for (j in i:n) {
# #       d <- mahalanobis(data[i, ], data[j, ], cov_matrix)
# #       dist_matrix[i, j] <- d
# #       dist_matrix[j, i] <- d
# #     }
# #   }
# #   
# #   # Convert to "dist" object
# #   return(as.dist(dist_matrix))
# # }
# # 
# # simulationFunction <-  function(tree, treeIdx){
# #   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
# #   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
# #   
# #   for(nC in numCharacters){
# #     for(sR in setRate){
# #       for(pc in proportionConflicting){
# #         for(d in numDatasets){
# #           
# #           datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #           pca <- prcomp(datasetD)
# #           plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #           pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #           njTree <- nj(pcDistMat)
# #           unrootedTree <- unroot(tree)
# #           unrootedNJTree <- unroot(njTree)
# #           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
# #           names(resVec) <- NULL
# #           resMat[nrow(resMat) + 1, ] <- resVec
# #         }
# #       }
# #     }
# #     for(vR in 1:3){
# #       for(pc in proportionConflicting){
# #         if(vR == 1){
# #           # 1, 10
# #           for(d in numDatasets){
# #             
# #             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #             pca <- prcomp(datasetD)
# #             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #             pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #             njTree <- nj(pcDistMat)
# #             unrootedTree <- unroot(tree)
# #             unrootedNJTree <- unroot(njTree)
# #             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
# #             names(resVec) <- NULL
# #             resMat[nrow(resMat) + 1, ] <- resVec
# #           }
# #         }else if (vR == 2){
# #           # 1, 1
# #           for(d in numDatasets){
# #             
# #             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #             pca <- prcomp(datasetD)
# #             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #             pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #             njTree <- nj(pcDistMat)
# #             unrootedTree <- unroot(tree)
# #             unrootedNJTree <- unroot(njTree)
# #             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
# #             names(resVec) <- NULL
# #             resMat[nrow(resMat) + 1, ] <- resVec
# #           }
# #         }else{
# #           # 10, 1
# #           for(d in numDatasets){
# #             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #             pca <- prcomp(datasetD)
# #             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #             pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #             njTree <- nj(pcDistMat)
# #             unrootedTree <- unroot(tree)
# #             unrootedNJTree <- unroot(njTree)
# #             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
# #             names(resVec) <- NULL
# #             resMat[nrow(resMat) + 1, ] <- resVec
# #           }
# #         }
# #       }
# #     }
# #   }
# #   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# # }
# # 
# # resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
# #   return(simulationFunction(phyloList[[i]], i))
# # })
# # 
# # saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistancePC1PC225char50char.rds")
# 
# 
# # simulation - 500 characters ---------------------------------------------
# 
# numCharacters <- c(500)
# 
# simulationFunction <- function(tree, treeIdx){
#   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
#   
#   for(nC in numCharacters){
#     for(sR in setRate){
#       for(pc in proportionConflicting){
#         for(d in numDatasets){
#           
#           datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#           rownames(datasetD) = tree$tip.label
#           for(t in 1:nC){
#             traits <- fastBM(tree, sig2 = sR)
#             datasetD[,t] <- traits
#           }
#           
#           whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#           for(t in whichConflicting){
#             newVec <- datasetD[,t]
#             tax <- sample(nrow(datasetD), 2, replace = F)
#             tax1Val <- newVec[tax[1]]
#             tax2Val <- newVec[tax[2]]
#             
#             datasetD[tax[1],t] <- tax2Val
#             datasetD[tax[2],t] <- tax1Val
#           }
#           dt <- as.data.table(datasetD)
#           rownames(dt) <- rownames(datasetD)
#           
#           data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#           pca <- prcomp(datasetD)
#           plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#           pcDistMat<- dist(plotdf)
#           njTree <- nj(pcDistMat)
#           unrootedTree <- unroot(tree)
#           unrootedNJTree <- unroot(njTree)
#           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
#           names(resVec) <- NULL
#           resMat[nrow(resMat) + 1, ] <- resVec
#         }
#       }
#     }
#     for(vR in 1:3){
#       for(pc in proportionConflicting){
#         if(vR == 1){
#           # 1, 10
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 1, scale = 10)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else if (vR == 2){
#           # 1, 1
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 1, scale = 1)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else{
#           # 10, 1
#           for(d in numDatasets){
#             
#             datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
#             rownames(datasetD) = tree$tip.label
#             for(t in 1:nC){
#               rate <- rgamma(1, shape = 10, scale = 1)
#               traits <- fastBM(tree, sig2 = rate)
#               datasetD[,t] <- traits
#             }
#             
#             whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
#             for(t in whichConflicting){
#               newVec <- datasetD[,t]
#               tax <- sample(nrow(datasetD), 2, replace = F)
#               tax1Val <- newVec[tax[1]]
#               tax2Val <- newVec[tax[2]]
#               
#               datasetD[tax[1],t] <- tax2Val
#               datasetD[tax[2],t] <- tax1Val
#             }
#             dt <- as.data.table(datasetD)
#             rownames(dt) <- rownames(datasetD)
#             
#             data.table::fwrite(dt,  paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = TRUE, verbose = FALSE)
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }
#       }
#     }
#   }
#   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# }
# 
# resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
#   return(simulationFunction(phyloList[[i]], i))
# })
# 
# saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults500char.rds")
# 
# simulationFunctionAllPCs <- function(tree, treeIdx){
#   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
#   
#   for(nC in numCharacters){
#     for(sR in setRate){
#       for(pc in proportionConflicting){
#         for(d in numDatasets){
#           datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#           pca <- prcomp(datasetD)
#           pcDistMat<- dist(pca$x)
#           njTree <- nj(pcDistMat)
#           unrootedTree <- unroot(tree)
#           unrootedNJTree <- unroot(njTree)
#           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
#           names(resVec) <- NULL
#           resMat[nrow(resMat) + 1, ] <- resVec
#         }
#       }
#     }
#     for(vR in 1:3){
#       for(pc in proportionConflicting){
#         if(vR == 1){
#           # 1, 10
#           for(d in numDatasets){
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else if (vR == 2){
#           # 1, 1
#           for(d in numDatasets){
#             
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }else{
#           # 10, 1
#           for(d in numDatasets){
#             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
#             pca <- prcomp(datasetD)
#             pcDistMat<- dist(pca$x)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
#             names(resVec) <- NULL
#             resMat[nrow(resMat) + 1, ] <- resVec
#           }
#         }
#       }
#     }
#   }
#   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# }
# 
# resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
#   return(simulationFunctionAllPCs(phyloList[[i]], i))
# })
# 
# saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsAllPCs500char.rds")
# 
# # mahalanobis_dist_matrix <- function(data) {
# #   # Ensure data is a matrix
# #   data <- as.matrix(data)
# #   
# #   # Compute the covariance matrix of the data
# #   cov_matrix <- cov(data)
# #   
# #   # Initialize distance matrix
# #   n <- nrow(data)
# #   dist_matrix <- matrix(0, n, n)
# #   rownames(dist_matrix) <- rownames(data)
# #   colnames(dist_matrix) <- rownames(data)
# #   
# #   # Compute pairwise Mahalanobis distances
# #   for (i in 1:n) {
# #     for (j in i:n) {
# #       d <- mahalanobis(data[i, ], data[j, ], cov_matrix)
# #       dist_matrix[i, j] <- d
# #       dist_matrix[j, i] <- d
# #     }
# #   }
# #   
# #   # Convert to "dist" object
# #   return(as.dist(dist_matrix))
# # }
# # 
# # simulationFunction <-  function(tree, treeIdx){
# #   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
# #   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
# #   
# #   for(nC in numCharacters){
# #     for(sR in setRate){
# #       for(pc in proportionConflicting){
# #         for(d in numDatasets){
# #           
# #           datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #           pca <- prcomp(datasetD)
# #           plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #           pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #           njTree <- nj(pcDistMat)
# #           unrootedTree <- unroot(tree)
# #           unrootedNJTree <- unroot(njTree)
# #           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
# #           names(resVec) <- NULL
# #           resMat[nrow(resMat) + 1, ] <- resVec
# #         }
# #       }
# #     }
# #     for(vR in 1:3){
# #       for(pc in proportionConflicting){
# #         if(vR == 1){
# #           # 1, 10
# #           for(d in numDatasets){
# #             
# #             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #             pca <- prcomp(datasetD)
# #             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #             pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #             njTree <- nj(pcDistMat)
# #             unrootedTree <- unroot(tree)
# #             unrootedNJTree <- unroot(njTree)
# #             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
# #             names(resVec) <- NULL
# #             resMat[nrow(resMat) + 1, ] <- resVec
# #           }
# #         }else if (vR == 2){
# #           # 1, 1
# #           for(d in numDatasets){
# #             
# #             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #             pca <- prcomp(datasetD)
# #             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #             pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #             njTree <- nj(pcDistMat)
# #             unrootedTree <- unroot(tree)
# #             unrootedNJTree <- unroot(njTree)
# #             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
# #             names(resVec) <- NULL
# #             resMat[nrow(resMat) + 1, ] <- resVec
# #           }
# #         }else{
# #           # 10, 1
# #           for(d in numDatasets){
# #             datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
# #             pca <- prcomp(datasetD)
# #             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
# #             pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
# #             njTree <- nj(pcDistMat)
# #             unrootedTree <- unroot(tree)
# #             unrootedNJTree <- unroot(njTree)
# #             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
# #             names(resVec) <- NULL
# #             resMat[nrow(resMat) + 1, ] <- resVec
# #           }
# #         }
# #       }
# #     }
# #   }
# #   return(list("phylo" = write.tree(tree), "resMat" = resMat))
# # }
# # 
# # resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
# #   return(simulationFunction(phyloList[[i]], i))
# # })
# # 
# # saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistancePC1PC2500char.rds")
# 

# LDDMM results -----------------------------------------------------------


### PC1 and 2 only ###

### non procrustes
#rSPRFunc
analysisFunction <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn

  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  pca <- prcomp(pcaMat)

  pcDistMat<- dist(pca$x[,1:2])
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistances/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), ".csv", sep = ""))
  
  return(resMat)
}

lddmmFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes",full.names = T)
doneFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistances")

#prune already analyzed lddmm files
pruneFunc <- function(lddmmFile){
  tag <- regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE))
  tag <- paste(tag, ".csv", sep = "")
  
  if(tag %in% doneFiles){
    
  }else{
    return(lddmmFile)
  }
  
}

lddmmFilesPruned <- parallel::mclapply(lddmmFiles, pruneFunc, mc.cores = 4)
lddmmFilesPruned <- unlist(lddmmFilesPruned)

rm(lddmmFiles)
rm(doneFiles)
gc()

res <- pbapply::pblapply(lddmmFilesPruned, FUN = analysisFunction)

resmat <- matrix(data = NA, nrow = length(res), ncol = 6)
for(i in 1:length(res)){
  resmat[i, ] <- as.numeric(res[[i]])
}
colnames(resmat) <- names(res[[i]])


### Procrustes
#rSPRFunc3
analysisFunctionProc <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  ### gpagen expects an array containing n x d slices where n is num landmarks and d is dimension of ladnmarks
  
  pcaMat <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  proc <- geomorph::gpagen(A = pcaMat, verbose = FALSE,print.progress = FALSE)
  
  pca <- geomorph::gm.prcomp(proc$coords)
  
  pcDistMat<- dist(pca$x[,1:2])
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc3(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustes/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), "Procrustes", ".csv", sep = ""))
  
  return(resMat)
}

lddmmFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes",full.names = T)
doneFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustes")

#prune already analyzed lddmm files
pruneFunc <- function(lddmmFile){
  tag <- regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE))
  tag <- paste(tag, "Procrustes.csv", sep = "")
  
  if(tag %in% doneFiles){
    
  }else{
    return(lddmmFile)
  }
  
}

lddmmFilesPruned <- parallel::mclapply(lddmmFiles, pruneFunc, mc.cores = 4)
lddmmFilesPruned <- unlist(lddmmFilesPruned)

rm(lddmmFiles)
rm(doneFiles)
gc()

res <- pbapply::pblapply(lddmmFilesPruned, FUN = analysisFunctionProc)

resmat <- matrix(data = NA, nrow = length(res), ncol = 6)
for(i in 1:length(res)){
  resmat[i, ] <- as.numeric(res[[i]])
}
colnames(resmat) <- names(res[[i]])

### analysis
doneFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustes/", full.names = T)

res <- parallel::mclapply(doneFiles, data.table::fread, mc.cores = 4)

resmat <- matrix(data = NA, nrow = length(res), ncol = 6)
for(i in 1:length(res)){
  resmat[i, ] <- as.numeric(res[[i]])
  if(i %% 100 == 0 ){
    print(i)
  }
}
colnames(resmat) <- names(res[[i]])



### ALL PCs

#rSPRFunc4
### non procrustes
analysisFunctionAll <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  pca <- prcomp(pcaMat)
  
  pcDistMat<- dist(pca$x)
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc4(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesAllPCs/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), ".csv", sep = ""))
  
  return(resMat)
}

lddmmFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes",full.names = T)
doneFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesAllPCs/")

#prune already analyzed lddmm files
pruneFunc <- function(lddmmFile){
  tag <- regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE))
  tag <- paste(tag, ".csv", sep = "")
  
  if(tag %in% doneFiles){
    
  }else{
    return(lddmmFile)
  }
  
}

lddmmFilesPruned <- parallel::mclapply(lddmmFiles, pruneFunc, mc.cores = 4)
lddmmFilesPruned <- unlist(lddmmFilesPruned)

rm(lddmmFiles)
rm(doneFiles)
gc()

res <- pbapply::pblapply(lddmmFilesPruned, FUN = analysisFunctionAll)

### Procrustes

#rSPRFunc2
analysisFunctionProcAll <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  ### gpagen expects an array containing n x d slices where n is num landmarks and d is dimension of ladnmarks
  
  pcaMat <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  proc <- geomorph::gpagen(A = pcaMat, verbose = FALSE,print.progress = FALSE)
  
  pca <- geomorph::gm.prcomp(proc$coords)
  
  pcDistMat<- dist(pca$x)
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc2(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustesAllPCs/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), "Procrustes", ".csv", sep = ""))
  
  return(resMat)
}

lddmmFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes",full.names = T)
doneFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustesAllPCs")

#prune already analyzed lddmm files
pruneFunc <- function(lddmmFile){
  tag <- regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE))
  tag <- paste(tag, "Procrustes.csv", sep = "")
  
  if(tag %in% doneFiles){
    
  }else{
    return(lddmmFile)
  }
  
}

lddmmFilesPruned <- parallel::mclapply(lddmmFiles, pruneFunc, mc.cores = 4)
lddmmFilesPruned <- unlist(lddmmFilesPruned)

rm(lddmmFiles)
rm(doneFiles)
gc()

res <- pbapply::pblapply(lddmmFilesPruned, FUN = analysisFunctionProcAll)



## function to find place in lddmm forward sims
maxTree <- function(lddmmFile, lm){
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  if(numLandmarks == lm){
    tree <- as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))
    return(tree)
  }
}

res <- parallel::mclapply(lddmmFiles, function(i){
  return(maxTree(i, 50))
}, mc.cores = 12)

x <- unlist(res)
table(x)


### efficient LDDMM proc ####

## PC1, 2
#rSPRFunc3
analysisFunctionProc <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  ### gpagen expects an array containing n x d slices where n is num landmarks and d is dimension of ladnmarks
  
  pcaMat <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  proc <- geomorph::gpagen(A = pcaMat, verbose = FALSE,print.progress = FALSE)
  
  pca <- geomorph::gm.prcomp(proc$coords)
  
  pcDistMat<- dist(pca$x[,1:2])
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc3(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustes/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), "Procrustes", ".csv", sep = ""))
  
  return(resMat)
}
#rSPRFunc4
analysisFunctionProc2 <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  ### gpagen expects an array containing n x d slices where n is num landmarks and d is dimension of ladnmarks
  
  pcaMat <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  proc <- geomorph::gpagen(A = pcaMat, verbose = FALSE,print.progress = FALSE)
  
  pca <- geomorph::gm.prcomp(proc$coords)
  
  pcDistMat<- dist(pca$x[,1:2])
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc4(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustes/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), "Procrustes", ".csv", sep = ""))
  
  return(resMat)
}

lddmmFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes",full.names = T)

doneFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustes")

#prune already analyzed lddmm files
pruneFunc <- function(lddmmFile){
  tag <- regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE))
  tag <- paste(tag, "Procrustes.csv", sep = "")
  
  if(tag %in% doneFiles){
    
  }else{
    return(lddmmFile)
  }
  
}

lddmmFilesPruned <- parallel::mclapply(lddmmFiles, pruneFunc, mc.cores = 8)
#lddmmFilesPruned <- unlist(lddmmFilesPruned)

rm(lddmmFiles)
rm(doneFiles)
gc()

data.table::fwrite(lddmmFilesPruned, "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/PC1PC2.tsv",col.names = FALSE)

rm(lddmmFilesPruned)
lddmmFilesPruned <- data.table::fread("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/PC1PC2.tsv", header = FALSE)
lddmmFilesPruned <- unlist(lddmmFilesPruned)

cutoff <-round(length(lddmmFilesPruned)/2)
lddmmFilesPrunedHalf1 <- lddmmFilesPruned[1:cutoff]
lddmmFilesPrunedHalf2 <- lddmmFilesPruned[(cutoff + 1) : length(lddmmFilesPruned)]

res <- pbapply::pblapply(lddmmFilesPrunedHalf1, FUN = analysisFunctionProc)
res <- pbapply::pblapply(lddmmFilesPrunedHalf2, FUN = analysisFunctionProc2)

# PCAll
#rSPRFunc2
analysisFunctionProcAll <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  ### gpagen expects an array containing n x d slices where n is num landmarks and d is dimension of ladnmarks
  
  pcaMat <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  proc <- geomorph::gpagen(A = pcaMat, verbose = FALSE,print.progress = FALSE)
  
  pca <- geomorph::gm.prcomp(proc$coords)
  
  pcDistMat<- dist(pca$x)
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc2(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustesAllPCs/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), "Procrustes", ".csv", sep = ""))
  
  return(resMat)
}
#rSPRFunc
analysisFunctionProcAll2 <- function(lddmmFile){
  tree <- phyloList[[1+as.numeric(regmatches(lddmmFile, regexpr("(?<=DimensionSimTreeIndex).*?(?=LM)", lddmmFile, perl = TRUE)))]]
  dimension <- (regmatches(lddmmFile, regexpr("(?<=LDDMMSimRes/).*?(?=Dimension)", lddmmFile, perl = TRUE)))
  numLandmarks <- as.numeric(regmatches(lddmmFile, regexpr("(?<=LM).*?(?=Alpha)", lddmmFile, perl = TRUE)))
  alpha <- as.numeric(regmatches(lddmmFile, regexpr("(?<=Alpha).*?(?=Dataset)", lddmmFile, perl = TRUE)))
  if(dimension == "three"){
    dimension <- 3
  }else if(dimension == "two"){
    dimension <- 2
  }
  
  landmarks <- read.delim(lddmmFile, header = F)
  cn <- colnames(landmarks)
  cn[1] <- "label"
  cn[2] <- "lm"
  colnames(landmarks) <- cn
  
  pcaMat <- matrix(data = NA, nrow = length(unique(landmarks$label)), ncol = dimension * numLandmarks)
  rownames(pcaMat) = unique(landmarks$label)
  for(i in unique(landmarks$label)){
    lm <- as.matrix(filter(landmarks, label == i))
    lm <- c(lm[, 3:ncol(lm)])
    for(j in 1:length(lm)){
      pcaMat[which(rownames(pcaMat)== i), j] = as.numeric(lm[j]) 
    }
  }
  ### gpagen expects an array containing n x d slices where n is num landmarks and d is dimension of ladnmarks
  
  pcaMat <- geomorph::arrayspecs(pcaMat, numLandmarks, k = dimension)
  proc <- geomorph::gpagen(A = pcaMat, verbose = FALSE,print.progress = FALSE)
  
  pca <- geomorph::gm.prcomp(proc$coords)
  
  pcDistMat<- dist(pca$x)
  njTree <- nj(pcDistMat)
  unrootedTree <- unroot(tree)
  unrootedNJTree <- unroot(njTree)
  
  ### this is the expensive step
  sprd <- rSPRFunc(unrootedNJTree, unrootedTree)
  
  resMat = as.data.frame(matrix(data = NA, nrow = 1, ncol = 6))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numLandmarks", "Dimension", "alpha")
  resMat[1,] <- c(dist.topo(unrootedNJTree, unrootedTree), sprd, sprd / Ntip(mast(unrootedNJTree, unrootedTree)),numLandmarks, dimension, alpha)
  resMat <- as.data.table(resMat)
  fwrite(resMat, paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustesAllPCs/", regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE)), "Procrustes", ".csv", sep = ""))
  
  return(resMat)
}

lddmmFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes",full.names = T)
doneFiles <- list.files("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMDistancesProcrustesAllPCs")

#prune already analyzed lddmm files
pruneFunc <- function(lddmmFile){
  tag <- regmatches(lddmmFile, regexpr("(?<=SimRes/LDDMMSimRes/).*?(?=.tsv)", lddmmFile, perl = TRUE))
  tag <- paste(tag, "Procrustes.csv", sep = "")
  
  if(tag %in% doneFiles){
    
  }else{
    return(lddmmFile)
  }
  
}

lddmmFilesPruned <- parallel::mclapply(lddmmFiles, pruneFunc, mc.cores = 4)
#lddmmFilesPruned <- unlist(lddmmFilesPruned)

rm(lddmmFiles)
rm(doneFiles)
gc()

data.table::fwrite(lddmmFilesPruned, "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/PCAll.tsv",col.names = FALSE)

rm(lddmmFilesPruned)
lddmmFilesPruned <- data.table::fread("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/PCAll.tsv", header = FALSE)
lddmmFilesPruned <- unlist(lddmmFilesPruned)

cutoff <-round(length(lddmmFilesPruned)/2)
lddmmFilesPrunedHalf1 <- lddmmFilesPruned[1:cutoff]
lddmmFilesPrunedHalf2 <- lddmmFilesPruned[(cutoff + 1) : length(lddmmFilesPruned)]

res <- pbapply::pblapply(lddmmFilesPrunedHalf1, FUN = analysisFunctionProcAll)
res <- pbapply::pblapply(lddmmFilesPrunedHalf2, FUN = analysisFunctionProcAll2)


# # Mongle et al. (2023) dataset --------------------------------------------
# 
# mea23 <- t(list2DF(read.nexus.data("/Users/levir/Documents/GitHub/PCAPhylogenetics/data/mongle_2023_rb.nex")))
# pcaMat <- matrix(data = NA, nrow = nrow(mea23), ncol = ncol(mea23))
# for(i in 1:ncol(mea23)){
#   pcaMat[,i] <- as.numeric(mea23[,i])
# }
# rownames(pcaMat) <- rownames(mea23)
# mea23PCA <- pcaMethods::pca(pcaMat, nPcs = 2, method = "ppca")
# 
# pcDistMat<- dist(mea23PCA@scores)
# njTree <- nj(pcDistMat)
# unrootedNJTree <- unroot(njTree)
# 
# treeSubset <- trees[sample(1:nrow(trees), 100000, replace = FALSE), ]
# phyloList2 <- list()
# for(i in 1:nrow(treeSubset)){
#   phyloList2[[i]] <- ape::read.tree(text = as.character(treeSubset[i, 5]))
# }
# 
# rfVec <- c()
# sprVec <- c()
# spreVec <- c()
# for(i in phyloList2){
#   tree <- i
#   rfVec <- c(rfVec, dist.topo(unrootedNJTree, unroot(tree)))
#   sprVec <- c(sprVec, rSPRFunc(unrootedNJTree, unroot(tree)))
#   spreVec <- c(spreVec, sprMast(unrootedNJTree, unroot(tree)))
# }
# print(sum(rfVec == 0))
# table(rfVec)
# 
# # Berger et al. (2015) dataset --------------------------------------------
# 
# bea2015 <- t(readxl::read_excel("/Users/levi/Documents/GitHub/PCAPhylogenetics/data/Berger_et_al_2015.xlsx"))
# bea2015 <- bea2015[2:13,]
# 
# ape::write.nexus.data(bea2015, "/Users/levi/Documents/GitHub/PCAPhylogenetics/data/Berger_et_al_2015.nex", format = "continuous")
# 
# pcaMat <- matrix(data = NA, nrow = nrow(bea2015), ncol = ncol(bea2015))
# for(i in 1:ncol(bea2015)){
#   pcaMat[,i] <- as.numeric(bea2015[,i])
# }
# bea15PCA <- pcaMethods::pca(pcaMat, nPcs = 2, method = "ppca")
# 
# plotdf <- data.frame(tax = rownames(bea2015), pc1 = bea15PCA@scores[,1], pc2 = bea15PCA@scores[,2])
# ggplot()+
#   geom_point(data = plotdf, aes(x = pc1, y = pc2, color = tax))+
#   geom_text(data = plotdf, aes(label = tax, x = pc1, y= pc2), hjust = 0, vjust = -0.5, size = 3) +
#   theme_minimal()
# 


# mahalanobis distance ----------------------------------------------------

mahalanobis_dist_matrix <- function(data) {
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Compute the covariance matrix of the data
  cov_matrix <- cov(data)
  
  # Initialize distance matrix
  n <- nrow(data)
  dist_matrix <- matrix(0, n, n)
  rownames(dist_matrix) <- rownames(data)
  colnames(dist_matrix) <- rownames(data)
  
  # Compute pairwise Mahalanobis distances
  for (i in 1:n) {
    for (j in i:n) {
      d <- mahalanobis(data[i, ], data[j, ], cov_matrix)
      dist_matrix[i, j] <- d
      dist_matrix[j, i] <- d
    }
  }
  
  # Convert to "dist" object
  return(as.dist(dist_matrix))
}

simulationFunction <-  function(tree, treeIdx){
  resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
  
  for(nC in numCharacters){
    for(sR in setRate){
      for(pc in proportionConflicting){
        for(d in numDatasets){
          
          datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "SetRate", sR,"PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
          pca <- prcomp(datasetD)
          plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
          pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
          njTree <- nj(pcDistMat)
          unrootedTree <- unroot(tree)
          unrootedNJTree <- unroot(njTree)
          resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
          names(resVec) <- NULL
          resMat[nrow(resMat) + 1, ] <- resVec
        }
      }
    }
    for(vR in 1:3){
      for(pc in proportionConflicting){
        if(vR == 1){
          # 1, 10
          for(d in numDatasets){
            
            datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else if (vR == 2){
          # 1, 1
          for(d in numDatasets){
            
            datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape1Scale1PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else{
          # 10, 1
          for(d in numDatasets){
            datasetD <- read.csv(paste("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/fastBMSimRes/fastBMTree", treeIdx, "NumChar", nC, "VarRatesShape10Scale10PropConflicting",pc,"Dataset",d,  ".csv", sep = ""), row.names = 1)
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), rSPRFunc(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }
      }
    }
  }
  return(list("phylo" = write.tree(tree), "resMat" = resMat))
}

resList <- pbapply::pblapply(1:length(phyloList), FUN = function(i){
  return(simulationFunction(phyloList[[i]], i))
})

saveRDS(resList, file = "/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistancePC1PC2.rds")

total <- 0
total0 <- 0
for(i in resList){
  df <- i[["resMat"]]
  total <- nrow(df) + total
  total0 <- sum(df$RF == 0) + total0  
}

rfDists <- c()
sprDists <- c()
for(i in resList){
  df <- i[["resMat"]]
  rfDists <- c(rfDists, df$RF)
  sprDists <- c(sprDists, df$SPR)
}
median(rfDists)
mean(sprDists)


# Brownian motion Analysis ----------------------------------------------------------------

simResPC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults.rds")
simRes2550PC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults25char50char.rds")
simRes500PC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults500char.rds")


simResPC1PC2 <- c(simResPC1PC2, simRes2550PC1PC2, simRes500PC1PC2)

for(i in 1:length(simResPC1PC2)){
  if(i == 1){
    concatDFPC1PC2 <- simResPC1PC2[[1]]$resMat
  }else{
    concatDFPC1PC2 <- rbind(concatDFPC1PC2, simResPC1PC2[[i]]$resMat)
  }
  if(i %% 10 == 0){
    print(i)
  }
}

concatDFPC1PC2$numCharacters <- as.factor(concatDFPC1PC2$numCharacters)
concatDFPC1PC2$setRate <- as.factor(concatDFPC1PC2$setRate)
concatDFPC1PC2$propConflicting <- as.factor(concatDFPC1PC2$propConflicting)
concatDFPC1PC2$varRateExpectation <- concatDFPC1PC2$variableRateShape / concatDFPC1PC2$variableRateScale
concatDFPC1PC2$varRateExpectation <- as.factor(concatDFPC1PC2$varRateExpectation)

#num trees identical
sum(concatDFPC1PC2$RF == 0) / nrow(concatDFPC1PC2)
sum(concatDFPC1PC2$RF == 0)
nrow(concatDFPC1PC2)
summary(concatDFPC1PC2$RF)
summary(concatDFPC1PC2$SPR)

for(i in unique(concatDFPC1PC2$numCharacters)){
  print(paste("nc: ", i, " with ", sum(filter(concatDFPC1PC2, numCharacters==i)$RF==0)))
}


simResAll <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsAllPCs.rds")
simRes2550All <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsAllPCs25char50char.rds")
simRes500All <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsAllPCs500char.rds")

simResAll <- c(simResAll, simRes2550All, simRes500All)

for(i in 1:length(simResAll)){
  if(i == 1){
    concatDFAll <- simResAll[[1]]$resMat
  }else{
    concatDFAll <- rbind(concatDFAll, simResAll[[i]]$resMat)
  }
  if(i %% 100 == 0){
    print(i)
  }
}

concatDFAll$numCharacters <- as.factor(concatDFAll$numCharacters)
concatDFAll$setRate <- as.factor(concatDFAll$setRate)
concatDFAll$propConflicting <- as.factor(concatDFAll$propConflicting)
concatDFAll$varRateExpectation <- concatDFAll$variableRateShape / concatDFAll$variableRateScale
concatDFAll$varRateExpectation <- as.factor(concatDFAll$varRateExpectation)

#num trees identical
sum(concatDFAll$RF == 0) / nrow(concatDFAll)
sum(concatDFAll$RF == 0) 
summary(concatDFAll$RF)
summary(concatDFAll$SPR)

for(i in unique(concatDFAll$numCharacters)){
  print(paste("nc: ", i, " with ", sum(filter(concatDFAll, numCharacters==i)$RF==0)))
}


simResMahalPC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistancePC1PC2.rds")
simResMahalPC1PC22550 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistancePC1PC225char50char.rds")
simResMahalPC1PC2500 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistancePC1PC2500char.rds")

simResMahalPC1PC2 <- c(simResMahalPC1PC2, simResMahalPC1PC22550, simResMahalPC1PC2500)

for(i in 1:length(simResMahalPC1PC2)){
  if(i == 1){
    concatDFMPC12 <- simResMahalPC1PC2[[1]]$resMat
  }else{
    concatDFMPC12 <- rbind(concatDFMPC12, simResMahalPC1PC2[[i]]$resMat)
  }
}
concatDFMPC12$numCharacters <- as.factor(concatDFMPC12$numCharacters)
concatDFMPC12$setRate <- as.factor(concatDFMPC12$setRate)
concatDFMPC12$propConflicting <- as.factor(concatDFMPC12$propConflicting)
concatDFMPC12$varRateExpectation <- concatDFMPC12$variableRateShape / concatDFMPC12$variableRateScale
concatDFMPC12$varRateExpectation <- as.factor(concatDFMPC12$varRateExpectation)

#num trees identical
sum(concatDFMPC12$RF == 0) / nrow(concatDFMPC12)
sum(concatDFMPC12$RF == 0)

summary(concatDFMPC12$SPR)
summary(concatDFMPC12$RF)


#are the distributions sig. different between set rate and varrate expectation across numcharacters
#spr
nc <- unique(concatDFPC1PC2$numCharacters)
varRatesLG <- TRUE
rate <- c(0.1, 1, 10)
meanDiff <- c()
for(c in nc){
  for(r in rate){
    print(paste(c,"characters for rate", r))
    print("Mean for set rate: ")
    print(mean(filter(concatDFPC1PC2, numCharacters == c, setRate == r, propConflicting ==0)$SPR))
    print("Mean for var rate: ")
    print(mean(filter(concatDFPC1PC2, numCharacters == c, varRateExpectation == r, propConflicting ==0)$SPR))
    
    if(mean(filter(concatDFPC1PC2, numCharacters == c, varRateExpectation == r, propConflicting ==0)$SPR) < mean(filter(concatDFPC1PC2, numCharacters == c, setRate == r, propConflicting ==0)$SPR)){
      varRatesLG <- FALSE
    }
    
    print(paste("p =", wilcox.test(filter(concatDFPC1PC2, numCharacters == c, setRate == r, propConflicting ==0)$SPR, filter(concatDFPC1PC2, numCharacters == c, varRateExpectation == r, propConflicting ==0)$SPR)$p.value))
    print("------------------------------")
  }
}
#rf
nc <- unique(concatDFPC1PC2$numCharacters)
varRatesLG <- TRUE
rate <- c(0.1, 1, 10)
meanDiff <- c()
for(c in nc){
  for(r in rate){
    print(paste(c,"characters for rate", r))
    print("Mean for set rate: ")
    print(mean(filter(concatDFPC1PC2, numCharacters == c, setRate == r, propConflicting ==0)$RF))
    print("Mean for var rate: ")
    print(mean(filter(concatDFPC1PC2, numCharacters == c, varRateExpectation == r, propConflicting ==0)$RF))
    
    if(mean(filter(concatDFPC1PC2, numCharacters == c, varRateExpectation == r, propConflicting ==0)$RF) < mean(filter(concatDFPC1PC2, numCharacters == c, setRate == r, propConflicting ==0)$RF)){
      varRatesLG <- FALSE
    }
    
    print(paste("p =", wilcox.test(filter(concatDFPC1PC2, numCharacters == c, setRate == r, propConflicting ==0)$SPR, filter(concatDFPC1PC2, numCharacters == c, varRateExpectation == r, propConflicting ==0)$SPR)$p.value))
    print("------------------------------")
  }
}

#spr
nc <- unique(concatDFAll$numCharacters)
rate <- c(0.1, 1, 10)
varRatesLG <- TRUE
for(c in nc){
  for(r in rate){
    print(paste(c,"characters for rate", r))
    print("Mean for set rate: ")
    print(mean(filter(concatDFAll, numCharacters == c, setRate == r, propConflicting ==0)$SPR))
    print("Mean for var rate: ")
    print(mean(filter(concatDFAll, numCharacters == c, varRateExpectation == r, propConflicting ==0)$SPR))
    if(mean(filter(concatDFAll, numCharacters == c, varRateExpectation == r, propConflicting ==0)$SPR) < mean(filter(concatDFAll, numCharacters == c, setRate == r, propConflicting ==0)$SPR)){
      varRatesLG <- FALSE
    }
    
    print(paste("p =", wilcox.test(filter(concatDFAll, numCharacters == c, setRate == r, propConflicting ==0)$SPR, filter(concatDFAll, numCharacters == c, varRateExpectation == r, propConflicting ==0)$SPR)$p.value))
    print("------------------------------")
  }
}
#RF
nc <- unique(concatDFAll$numCharacters)
rate <- c(0.1, 1, 10)
varRatesLG <- TRUE
for(c in nc){
  for(r in rate){
    print(paste(c,"characters for rate", r))
    print("Mean for set rate: ")
    print(mean(filter(concatDFAll, numCharacters == c, setRate == r, propConflicting ==0)$RF))
    print("Mean for var rate: ")
    print(mean(filter(concatDFAll, numCharacters == c, varRateExpectation == r, propConflicting ==0)$RF))
    if(mean(filter(concatDFAll, numCharacters == c, varRateExpectation == r, propConflicting ==0)$RF) < mean(filter(concatDFAll, numCharacters == c, setRate == r, propConflicting ==0)$RF)){
      varRatesLG <- FALSE
    }
    
    print(paste("p =", wilcox.test(filter(concatDFAll, numCharacters == c, setRate == r, propConflicting ==0)$RF, filter(concatDFAll, numCharacters == c, varRateExpectation == r, propConflicting ==0)$RF)$p.value))
    print("------------------------------")
  }
}



nc <- unique(concatDFPC1PC2$numCharacters)
pc <- unique(concatDFPC1PC2$propConflicting)
meanDiff <- c()
for(c in nc){
  for(p in 2:length(pc)){
    print(paste(c,"characters"))
    print(paste("Mean for pc:", pc[p]))
    print(mean(filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p])$SPR))
    print(paste("Mean for pc: ", pc[p-1]))
    print(mean(filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p-1])$SPR))
    
    print(paste("p =", wilcox.test(filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p])$SPR, filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p-1])$SPR)$p.value))
    print("------------------------------")
  }
}

nc <- unique(concatDFAll$numCharacters)
pc <- unique(concatDFAll$propConflicting)
meanDiff <- c()
for(c in nc){
  for(p in 2:length(pc)){
    print(paste(c,"characters"))
    print(paste("Mean for pc:", pc[p]))
    print(mean(filter(concatDFAll, numCharacters == c, propConflicting ==pc[p])$SPR))
    print(paste("Mean for pc: ", pc[p-1]))
    print(mean(filter(concatDFAll, numCharacters == c, propConflicting ==pc[p-1])$SPR))
    
    print(paste("p =", wilcox.test(filter(concatDFAll, numCharacters == c, propConflicting ==pc[p])$SPR, filter(concatDFAll, numCharacters == c, propConflicting ==pc[p-1])$SPR)$p.value))
    print("------------------------------")
  }
}

nc <- unique(concatDFPC1PC2$numCharacters)
pc <- unique(concatDFPC1PC2$propConflicting)
meanDiff <- c()
for(c in nc){
  for(p in 2:length(pc)){
    print(paste(c,"characters"))
    print(paste("Mean for pc:", pc[p]))
    print(mean(filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p])$RF))
    print(paste("Mean for pc: ", pc[p-1]))
    print(mean(filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p-1])$RF))
    
    print(paste("p =", wilcox.test(filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p])$RF, filter(concatDFPC1PC2, numCharacters == c, propConflicting ==pc[p-1])$RF)$p.value))
    print("------------------------------")
  }
}

nc <- unique(concatDFAll$numCharacters)
pc <- unique(concatDFAll$propConflicting)
meanDiff <- c()
for(c in nc){
  for(p in 2:length(pc)){
    print(paste(c,"characters"))
    print(paste("Mean for pc:", pc[p]))
    print(mean(filter(concatDFAll, numCharacters == c, propConflicting ==pc[p])$RF))
    print(paste("Mean for pc: ", pc[p-1]))
    print(mean(filter(concatDFAll, numCharacters == c, propConflicting ==pc[p-1])$RF))
    
    print(paste("p =", wilcox.test(filter(concatDFAll, numCharacters == c, propConflicting ==pc[p])$RF, filter(concatDFAll, numCharacters == c, propConflicting ==pc[p-1])$RF)$p.value))
    print("------------------------------")
  }
}

nc <- unique(concatDFPC1PC2$numCharacters)
for(c in nc){
  print(paste("Num trees identical for char", c))
  print(nrow(filter(concatDFPC1PC2, RF==0, numCharacters == c)))
}

nc <- unique(concatDFAll$numCharacters)
for(c in nc){
  print(paste("Num trees identical for char", c))
  print(nrow(filter(concatDFAll, RF==0, numCharacters == c)))
}

# LDDMM analysis ----------------------------------------------------------

resPC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/lddmmPC12results.rds")
resAll <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/lddmmPCallresults.rds")
resProcPC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/lddmmProcrustesPC12results.rds")
resProcAll <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/lddmmProcrustesPCallresults.rds")

sum(resPC1PC2$RF == 0) / nrow(resPC1PC2)
sum(resPC1PC2$RF == 0)
summary(resPC1PC2$RF)
summary(resPC1PC2$SPR)

sum(resAll$RF == 0) / nrow(resAll)
sum(resAll$RF == 0)
summary(resAll$RF)
summary(resAll$SPR)

wilcox.test(filter(resPC1PC2, numLandmarks == 10)$SPR, filter(resPC1PC2, numLandmarks == 25)$SPR)
mean(filter(resPC1PC2, numLandmarks == 10)$SPR)
mean(filter(resPC1PC2, numLandmarks == 25)$SPR)

wilcox.test(filter(resPC1PC2, numLandmarks == 25)$SPR, filter(resPC1PC2, numLandmarks == 50)$SPR)
mean(filter(resPC1PC2, numLandmarks == 25)$SPR)
mean(filter(resPC1PC2, numLandmarks == 50)$SPR)

wilcox.test(filter(resPC1PC2, numLandmarks == 10)$SPR, filter(resPC1PC2, numLandmarks == 50)$SPR)
mean(filter(resPC1PC2, numLandmarks == 10)$SPR)
mean(filter(resPC1PC2, numLandmarks == 50)$SPR)

wilcox.test(filter(resAll, numLandmarks == 10)$SPR, filter(resAll, numLandmarks == 25)$SPR)
mean(filter(resAll, numLandmarks == 10)$SPR)
mean(filter(resAll, numLandmarks == 25)$SPR)

wilcox.test(filter(resAll, numLandmarks == 25)$SPR, filter(resAll, numLandmarks == 50)$SPR)
mean(filter(resAll, numLandmarks == 25)$SPR)
mean(filter(resAll, numLandmarks == 50)$SPR)

wilcox.test(filter(resAll, numLandmarks == 10)$SPR, filter(resAll, numLandmarks == 50)$SPR)
mean(filter(resAll, numLandmarks == 10)$SPR)
mean(filter(resAll, numLandmarks == 50)$SPR)


wilcox.test(filter(resPC1PC2, numLandmarks == 10)$RF, filter(resPC1PC2, numLandmarks == 25)$RF)
mean(filter(resPC1PC2, numLandmarks == 10)$RF)
mean(filter(resPC1PC2, numLandmarks == 25)$RF)

wilcox.test(filter(resPC1PC2, numLandmarks == 25)$RF, filter(resPC1PC2, numLandmarks == 50)$RF)
mean(filter(resPC1PC2, numLandmarks == 25)$RF)
mean(filter(resPC1PC2, numLandmarks == 50)$RF)

wilcox.test(filter(resPC1PC2, numLandmarks == 10)$RF, filter(resPC1PC2, numLandmarks == 50)$RF)
mean(filter(resPC1PC2, numLandmarks == 10)$RF)
mean(filter(resPC1PC2, numLandmarks == 50)$RF)

wilcox.test(filter(resAll, numLandmarks == 10)$RF, filter(resAll, numLandmarks == 25)$RF)
mean(filter(resAll, numLandmarks == 10)$RF)
mean(filter(resAll, numLandmarks == 25)$RF)

wilcox.test(filter(resAll, numLandmarks == 25)$RF, filter(resAll, numLandmarks == 50)$RF)
mean(filter(resAll, numLandmarks == 25)$RF)
mean(filter(resAll, numLandmarks == 50)$RF)

wilcox.test(filter(resAll, numLandmarks == 10)$RF, filter(resAll, numLandmarks == 50)$RF)
mean(filter(resAll, numLandmarks == 10)$RF)
mean(filter(resAll, numLandmarks == 50)$RF)

wilcox.test(filter(resPC1PC2, Dimension == 2)$SPR, filter(resPC1PC2, Dimension == 3)$SPR)
mean(filter(resPC1PC2, Dimension == 2)$SPR)
mean(filter(resPC1PC2, Dimension == 3)$SPR)

wilcox.test(filter(resAll, Dimension == 2)$SPR, filter(resAll, Dimension == 3)$SPR)
mean(filter(resAll, Dimension == 2)$SPR)
mean(filter(resAll, Dimension == 3)$SPR)

wilcox.test(filter(resPC1PC2, Dimension == 2)$RF, filter(resPC1PC2, Dimension == 3)$RF)
mean(filter(resPC1PC2, Dimension == 2)$RF)
mean(filter(resPC1PC2, Dimension == 3)$RF)

wilcox.test(filter(resAll, Dimension == 2)$RF, filter(resAll, Dimension == 3)$RF)
mean(filter(resAll, Dimension == 2)$RF)
mean(filter(resAll, Dimension == 3)$RF)


sum(resProcPC1PC2$RF == 0) / nrow(resProcPC1PC2)
sum(resProcPC1PC2$RF == 0)
summary(resProcPC1PC2$RF)
summary(resProcPC1PC2$SPR)

sum(resProcAll$RF == 0) / nrow(resProcAll)
sum(resProcAll$RF == 0)
summary(resProcAll$RF)
summary(resProcAll$SPR)

wilcox.test(filter(resProcPC1PC2, numLandmarks == 10)$SPR, filter(resProcPC1PC2, numLandmarks == 25)$SPR)
mean(filter(resProcPC1PC2, numLandmarks == 10)$SPR)
mean(filter(resProcPC1PC2, numLandmarks == 25)$SPR)

wilcox.test(filter(resProcPC1PC2, numLandmarks == 25)$SPR, filter(resProcPC1PC2, numLandmarks == 50)$SPR)
mean(filter(resProcPC1PC2, numLandmarks == 25)$SPR)
mean(filter(resProcPC1PC2, numLandmarks == 50)$SPR)

wilcox.test(filter(resProcPC1PC2, numLandmarks == 10)$SPR, filter(resProcPC1PC2, numLandmarks == 50)$SPR)
mean(filter(resProcPC1PC2, numLandmarks == 10)$SPR)
mean(filter(resProcPC1PC2, numLandmarks == 50)$SPR)

wilcox.test(filter(resProcAll, numLandmarks == 10)$SPR, filter(resProcAll, numLandmarks == 25)$SPR)
mean(filter(resProcAll, numLandmarks == 10)$SPR)
mean(filter(resProcAll, numLandmarks == 25)$SPR)

wilcox.test(filter(resProcAll, numLandmarks == 25)$SPR, filter(resProcAll, numLandmarks == 50)$SPR)
mean(filter(resProcAll, numLandmarks == 25)$SPR)
mean(filter(resProcAll, numLandmarks == 50)$SPR)

wilcox.test(filter(resProcAll, numLandmarks == 10)$SPR, filter(resProcAll, numLandmarks == 50)$SPR)
mean(filter(resProcAll, numLandmarks == 10)$SPR)
mean(filter(resProcAll, numLandmarks == 50)$SPR)


wilcox.test(filter(resProcPC1PC2, numLandmarks == 10)$RF, filter(resProcPC1PC2, numLandmarks == 25)$RF)
mean(filter(resProcPC1PC2, numLandmarks == 10)$RF)
mean(filter(resProcPC1PC2, numLandmarks == 25)$RF)

wilcox.test(filter(resProcPC1PC2, numLandmarks == 25)$RF, filter(resProcPC1PC2, numLandmarks == 50)$RF)
mean(filter(resProcPC1PC2, numLandmarks == 25)$RF)
mean(filter(resProcPC1PC2, numLandmarks == 50)$RF)

wilcox.test(filter(resProcPC1PC2, numLandmarks == 10)$RF, filter(resProcPC1PC2, numLandmarks == 50)$RF)
mean(filter(resProcPC1PC2, numLandmarks == 10)$RF)
mean(filter(resProcPC1PC2, numLandmarks == 50)$RF)

wilcox.test(filter(resProcAll, numLandmarks == 10)$RF, filter(resProcAll, numLandmarks == 25)$RF)
mean(filter(resProcAll, numLandmarks == 10)$RF)
mean(filter(resProcAll, numLandmarks == 25)$RF)

wilcox.test(filter(resProcAll, numLandmarks == 25)$RF, filter(resProcAll, numLandmarks == 50)$RF)
mean(filter(resProcAll, numLandmarks == 25)$RF)
mean(filter(resProcAll, numLandmarks == 50)$RF)

wilcox.test(filter(resProcAll, numLandmarks == 10)$RF, filter(resProcAll, numLandmarks == 50)$RF)
mean(filter(resProcAll, numLandmarks == 10)$RF)
mean(filter(resProcAll, numLandmarks == 50)$RF)

wilcox.test(filter(resProcPC1PC2, Dimension == 2)$SPR, filter(resProcPC1PC2, Dimension == 3)$SPR)
mean(filter(resProcPC1PC2, Dimension == 2)$SPR)
mean(filter(resProcPC1PC2, Dimension == 3)$SPR)

wilcox.test(filter(resProcAll, Dimension == 2)$SPR, filter(resProcAll, Dimension == 3)$SPR)
mean(filter(resProcAll, Dimension == 2)$SPR)
mean(filter(resProcAll, Dimension == 3)$SPR)

wilcox.test(filter(resProcPC1PC2, Dimension == 2)$RF, filter(resProcPC1PC2, Dimension == 3)$RF)
mean(filter(resProcPC1PC2, Dimension == 2)$RF)
mean(filter(resProcPC1PC2, Dimension == 3)$RF)

wilcox.test(filter(resProcAll, Dimension == 2)$RF, filter(resProcAll, Dimension == 3)$RF)
mean(filter(resProcAll, Dimension == 2)$RF)
mean(filter(resProcAll, Dimension == 3)$RF)

