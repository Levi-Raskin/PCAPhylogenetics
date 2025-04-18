# Setup -------------------------------------------------------------------

sprMast <- function(tree1, tree2){
  sprDist <- phangorn::SPR.dist(tree1, tree2)
  mastSize <- Ntip(phangorn::mast(tree1, tree2))
  return(sprDist/mastSize)
}

gatree <- read.newick(text = "(Papio:0,((Colobus:0.075065,Hylobates:0.178754):0.125803,(Pongo_pygmaeus:0.325375,(Gorilla_gorilla:0.144516,(Pan_troglodytes:0.083517,(Sahelanthropus_tchadensis:0.244099,(Ardipithecus_ramidus:0.091594,(Australopithecus_anamensis:0.212352,(Australopithecus_afarensis:0.108627,(((Kenyanthropus_platyops:0.112281,((Homo_habilis:0.204068,(Homo_ergaster:0.035543,Homo_sapiens:0.175628):0.090888):0.243237,Homo_rudolfensis:0.062298):0.135252):0.127238,Australopithecus_africanus:0.245198):0.266149,(Australopithecus_garhi:0.032466,(Paranthropus_aethiopicus:4.60512e-06,(Paranthropus_robustus:0.243695,Paranthropus_boisei:0.103549):0.296051):0.764347):0.038617):0.067820):0.159549):0.150527):0.091403):0.133450):0.091660):0.089983):0.213719):0.111894);")
plot(gatree)

# Parameters --------------------------------------------------------------

numSims <- 1000
numCharacters <- 100
numHomoplasticChars <- 10

# Simulation --------------------------------------------------------------

spre <- c()
rf <- c()
sprdist <- c()

for(j in 1:numSims){
  charmat <- matrix(data = NA, nrow = Ntip(gatree), ncol = numCharacters)
  rownames(charmat) = gatree$tip.label
  for(i in 1:numCharacters){
    traits <- fastBM(gatree, sig2 = 1)
    #traits <- as.numeric(simSeq(gatree, l = 1))
    charmat[,i] <- traits
  }
  
  pca <- prcomp(charmat)
  plotdf <- data.frame(tax = rownames(charmat), pc1 = pca$x[,1], pc2 = pca$x[,2])
  
  pcDistMat<- dist(plotdf)
  njTree <- nj(pcDistMat)
  #plot(njTree)
  spre <- c(spre, sprMast(unroot(gatree), njTree))
  sprdist <- c(sprdist, TreeDist::SPRDist(unroot(gatree), njTree))
  rf <- c(rf, dist.topo(unroot(gatree), njTree))
  print(paste("NG: ", j ," RF:" , dist.topo(unroot(gatree), njTree)))
}

ggplot()+
  geom_point(data = plotdf, aes(x = pc1, y = pc2, color = tax))+
  geom_text(data = plotdf, aes(label = tax, x = pc1, y= pc2), hjust = 0, vjust = -0.5, size = 3) +
  theme_minimal()

hist(sprdist)