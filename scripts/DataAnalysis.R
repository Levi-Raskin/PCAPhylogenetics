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
trees <- data.table::fread("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/hominin.trees")
#burn in
trees <- trees[round(0.1 * nrow(trees)) : nrow(trees),]
# sample 1,000 trees from the posterior distribution
treeSubset <- trees[sample(1:nrow(trees), 1000, replace = FALSE), ]
#convert to phylo object
phyloList <- list()
for(i in 1:nrow(treeSubset)){
  phyloList[[i]] <- ape::read.tree(text = as.character(treeSubset[i, 5]))
}


# Parameters --------------------------------------------------------------

numDatasets <- 100
numCharacters <- c(10, 100, 250)
setRate <- c(0.1, 1, 10)
proportionConflicting <- c(0, 0.1, 0.25, 0.5) #proportion of characters that exhibit conflicting phylogenetic signal
numTaxaConflicting <- c(2) #just rearranging 2 taxa to start
#variableRates # drawn from gamma with parameters (1,10); (1,1); (10,1)

# Simulation --------------------------------------------------------------

### Non parallelized:

# resList <- list()
# 
# for(i in 1:length(phyloList)){
#   tree <- phyloList[[i]]
#   
#   #resmat:
#   #columns: rf, spr, spre, numCharacters, setRate, variableRateShape, variableRateScale
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
#           pca <- prcomp(datasetD)
#           plotdf <- data.frame(tax = rownames(datasetD), pc1 = pca$x[,1], pc2 = pca$x[,2])
#           pcDistMat<- dist(plotdf)
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
#               rate <- rgamma(1, shape = 1, rate = 10)
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
#               rate <- rgamma(1, shape = 1, rate = 1)
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
#               rate <- rgamma(1, shape = 10, rate = 1)
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
#   
#   resList[[i]] <- list("phylo" = write.tree(tree), "resMat" = resMat)
#   print(i)
# }


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
              rate <- rgamma(1, shape = 1, rate = 10)
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
              rate <- rgamma(1, shape = 1, rate = 1)
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
              rate <- rgamma(1, shape = 10, rate = 1)
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
median(sprDists)
### sample multiple individuals per taxon

# Mongle et al. (2023) dataset --------------------------------------------

mea23 <- t(list2DF(read.nexus.data("/Users/levi/Documents/GitHub/PCAPhylogenetics/data/mongle_2023_rb.nex")))
pcaMat <- matrix(data = NA, nrow = nrow(mea23), ncol = ncol(mea23))
for(i in 1:ncol(mea23)){
  pcaMat[,i] <- as.numeric(mea23[,i])
}
mea23PCA <- pcaMethods::pca(pcaMat, nPcs = 2, method = "ppca")

plotdf <- data.frame(tax = rownames(mea23), pc1 = mea23PCA@scores[,1], pc2 = mea23PCA@scores[,2])
ggplot()+
  geom_point(data = plotdf, aes(x = pc1, y = pc2, color = tax))+
  geom_text(data = plotdf, aes(label = tax, x = pc1, y= pc2), hjust = 0, vjust = -0.5, size = 3) +
  theme_minimal()

hist(sprdist)


# Berger et al. (2015) dataset --------------------------------------------

bea2015 <- t(readxl::read_excel("/Users/levi/Documents/GitHub/PCAPhylogenetics/data/Berger_et_al_2015.xlsx"))
bea2015 <- bea2015[2:13,]

ape::write.nexus.data(bea2015, "/Users/levi/Documents/GitHub/PCAPhylogenetics/data/Berger_et_al_2015.nex", format = "continuous")

pcaMat <- matrix(data = NA, nrow = nrow(bea2015), ncol = ncol(bea2015))
for(i in 1:ncol(bea2015)){
  pcaMat[,i] <- as.numeric(bea2015[,i])
}
bea15PCA <- pcaMethods::pca(pcaMat, nPcs = 2, method = "ppca")

plotdf <- data.frame(tax = rownames(bea2015), pc1 = bea15PCA@scores[,1], pc2 = bea15PCA@scores[,2])
ggplot()+
  geom_point(data = plotdf, aes(x = pc1, y = pc2, color = tax))+
  geom_text(data = plotdf, aes(label = tax, x = pc1, y= pc2), hjust = 0, vjust = -0.5, size = 3) +
  theme_minimal()

