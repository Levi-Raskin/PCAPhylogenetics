# Setup -------------------------------------------------------------------
library(ape) 
library(ggplot2)
library(phangorn)
library(phytools)

set.seed(5)

sprMast <- function(tree1, tree2){
  sprDist <- phangorn::SPR.dist(tree1, tree2)
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
treeSubset <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/sampledTrees.csv")

#convert to phylo object
phyloList <- list()
for(i in 1:nrow(treeSubset)){
  phyloList[[i]] <- ape::read.tree(text = as.character(treeSubset[i, 2]))
}


# Parameters --------------------------------------------------------------

numDatasets <- 100
numCharacters <- c(10, 100, 250)
setRate <- c(0.1, 1, 10)
proportionConflicting <- c(0, 0.1, 0.25, 0.5) #proportion of characters that exhibit conflicting phylogenetic signal
# proportionConflicting <- c(0)
numTaxaConflicting <- c(2) #just rearranging 2 taxa to start
#variableRates # drawn from gamma with parameters (1,10); (1,1); (10,1)

# Simulation --------------------------------------------------------------

simulationFunction <- function(tree){
  resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
  
  for(nC in numCharacters){
    for(sR in setRate){
      for(pc in proportionConflicting){
        for(d in numDatasets){
          
          datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
          rownames(datasetD) = tree$tip.label
          for(t in 1:nC){
            traits <- fastBM(tree, sig2 = sR)
            datasetD[,t] <- traits
          }
          
          whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
          for(t in whichConflicting){
            newVec <- datasetD[,t]
            tax <- sample(nrow(datasetD), 2, replace = F)
            tax1Val <- newVec[tax[1]]
            tax2Val <- newVec[tax[2]]
            
            datasetD[tax[1],t] <- tax2Val
            datasetD[tax[2],t] <- tax1Val
          }
          
          pca <- prcomp(datasetD)
          plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
          pcDistMat<- dist(plotdf)
          njTree <- nj(pcDistMat)
          unrootedTree <- unroot(tree)
          unrootedNJTree <- unroot(njTree)
          resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
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
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 10)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- dist(plotdf)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else if (vR == 2){
          # 1, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- dist(plotdf)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else{
          # 10, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 10, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- dist(plotdf)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }
      }
    }
  }
  return(list("phylo" = write.tree(tree), "resMat" = resMat))
}

resList <- pbapply::pblapply(phyloList, simulationFunction)

saveRDS(resList, file = "Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults.rds")

simulationFunctionAllPCs <- function(tree){
  resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
  
  for(nC in numCharacters){
    for(sR in setRate){
      for(pc in proportionConflicting){
        for(d in numDatasets){
          
          datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
          rownames(datasetD) = tree$tip.label
          for(t in 1:nC){
            traits <- fastBM(tree, sig2 = sR)
            datasetD[,t] <- traits
          }
          
          whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
          for(t in whichConflicting){
            newVec <- datasetD[,t]
            tax <- sample(nrow(datasetD), 2, replace = F)
            tax1Val <- newVec[tax[1]]
            tax2Val <- newVec[tax[2]]
            
            datasetD[tax[1],t] <- tax2Val
            datasetD[tax[2],t] <- tax1Val
          }
          
          pca <- prcomp(datasetD)
          pcDistMat<- dist(pca$x)
          njTree <- nj(pcDistMat)
          unrootedTree <- unroot(tree)
          unrootedNJTree <- unroot(njTree)
          resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
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
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 10)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            pcDistMat<- dist(pca$x)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else if (vR == 2){
          # 1, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            pcDistMat<- dist(pca$x)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else{
          # 10, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 10, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            pcDistMat<- dist(pca$x)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }
      }
    }
  }
  return(list("phylo" = write.tree(tree), "resMat" = resMat))
}

resList <- pbapply::pblapply(phyloList, simulationFunctionAllPCs)

saveRDS(resList, file = "Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsAllPCs.rds")

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

### sample multiple individuals per taxon
# simulationFunction <- function(tree){
#   resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
#   colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
#   
#   for(nC in numCharacters){
#     for(sR in setRate){
#       for(pc in proportionConflicting){
#         for(d in numDatasets){
#           
#           sim <- Rphylopars::simtraits(
#             ntraits      = nC,         # number of traits
#             nreps        = 10,        # individuals per taxon
#             tree         = tree,      # NULL → generates a random tree of size ntaxa
#             intraspecific = 0.1,     # within-species variance (per individual)
#             model        = "BM"       # Brownian motion
#           )
#           datasetD = as.matrix(sim$trait_data[,2:(nC+1)])
#           rownames(datasetD) = sim$trait_data[,1]
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
#           pca <- prcomp(datasetD)
#           pcDistMat<- dist(pca$x[,1:2])
#           njTree <- nj(pcDistMat)
#           unrootedTree <- unroot(tree)
#           unrootedNJTree <- unroot(njTree)
#           resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
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
#             
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
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
#             
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
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
#             
#             pca <- prcomp(datasetD)
#             plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#             pcDistMat<- dist(plotdf)
#             njTree <- nj(pcDistMat)
#             unrootedTree <- unroot(tree)
#             unrootedNJTree <- unroot(njTree)
#             resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
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
# resList <- pbapply::pblapply(phyloList, simulationFunction)
# 
# saveRDS(resList, file = "Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsIntraspecificVar.rds")
# 
# ggplot()+
#   geom_point(data= plotdf, aes(x=pc1, y = pc2, color = tax))
# 
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
#   sprVec <- c(sprVec, SPR.dist(unrootedNJTree, unroot(tree)))
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

simulationFunction <- function(tree){
  resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
  
  for(nC in numCharacters){
    for(sR in setRate){
      for(pc in proportionConflicting){
        for(d in numDatasets){
          
          datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
          rownames(datasetD) = tree$tip.label
          for(t in 1:nC){
            traits <- fastBM(tree, sig2 = sR)
            datasetD[,t] <- traits
          }
          
          whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
          for(t in whichConflicting){
            newVec <- datasetD[,t]
            tax <- sample(nrow(datasetD), 2, replace = F)
            tax1Val <- newVec[tax[1]]
            tax2Val <- newVec[tax[2]]
            
            datasetD[tax[1],t] <- tax2Val
            datasetD[tax[2],t] <- tax1Val
          }
          
          pca <- prcomp(datasetD)
          plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
          pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
          njTree <- nj(pcDistMat)
          unrootedTree <- unroot(tree)
          unrootedNJTree <- unroot(njTree)
          resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
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
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 10)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else if (vR == 2){
          # 1, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else{
          # 10, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 10, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
            pcDistMat<- mahalanobis_dist_matrix(plotdf[,2:3])
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }
      }
    }
  }
  return(list("phylo" = write.tree(tree), "resMat" = resMat))
}

resList <- pbapply::pblapply(phyloList, simulationFunction)

saveRDS(resList, file = "Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistancePC1PC2.rds")

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

simulationFunctionAllPCs <- function(tree){
  resMat <- data.frame(matrix(data = NA, nrow = 0, ncol = 8))
  colnames(resMat) <- c("RF", "SPR", "SPRe", "numCharacters", "setRate", "variableRateShape", "variableRateScale", "propConflicting")
  
  for(nC in numCharacters){
    for(sR in setRate){
      for(pc in proportionConflicting){
        for(d in numDatasets){
          
          datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
          rownames(datasetD) = tree$tip.label
          for(t in 1:nC){
            traits <- fastBM(tree, sig2 = sR)
            datasetD[,t] <- traits
          }
          
          whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
          for(t in whichConflicting){
            newVec <- datasetD[,t]
            tax <- sample(nrow(datasetD), 2, replace = F)
            tax1Val <- newVec[tax[1]]
            tax2Val <- newVec[tax[2]]
            
            datasetD[tax[1],t] <- tax2Val
            datasetD[tax[2],t] <- tax1Val
          }
          
          pca <- prcomp(datasetD)
          pcDistMat<- mahalanobis_dist_matrix(pca$x)
          njTree <- nj(pcDistMat)
          unrootedTree <- unroot(tree)
          unrootedNJTree <- unroot(njTree)
          resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, sR, NA, NA, pc) 
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
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 10)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            pcDistMat<- mahalanobis_dist_matrix(pca$x)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 10, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else if (vR == 2){
          # 1, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 1, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            pcDistMat<- mahalanobis_dist_matrix(pca$x)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 1, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }else{
          # 10, 1
          for(d in numDatasets){
            
            datasetD <- matrix(data = NA, nrow = Ntip(tree), ncol = nC)
            rownames(datasetD) = tree$tip.label
            for(t in 1:nC){
              rate <- rgamma(1, shape = 10, scale = 1)
              traits <- fastBM(tree, sig2 = rate)
              datasetD[,t] <- traits
            }
            
            whichConflicting <- sample(1:nC, size = round(pc * nC), replace = FALSE)
            for(t in whichConflicting){
              newVec <- datasetD[,t]
              tax <- sample(nrow(datasetD), 2, replace = F)
              tax1Val <- newVec[tax[1]]
              tax2Val <- newVec[tax[2]]
              
              datasetD[tax[1],t] <- tax2Val
              datasetD[tax[2],t] <- tax1Val
            }
            
            pca <- prcomp(datasetD)
            pcDistMat<- mahalanobis_dist_matrix(pca$x)
            njTree <- nj(pcDistMat)
            unrootedTree <- unroot(tree)
            unrootedNJTree <- unroot(njTree)
            resVec <- c(dist.topo(unrootedNJTree, unrootedTree), SPR.dist(unrootedNJTree, unrootedTree), sprMast(unrootedNJTree, unrootedTree), nC, NA, 10, 1, pc) 
            names(resVec) <- NULL
            resMat[nrow(resMat) + 1, ] <- resVec
          }
        }
      }
    }
  }
  return(list("phylo" = write.tree(tree), "resMat" = resMat))
}

resList <- pbapply::pblapply(phyloList, simulationFunctionAllPCs)

saveRDS(resList, file = "Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResultsMahalanobisDistanceAllPCs.rds")



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
