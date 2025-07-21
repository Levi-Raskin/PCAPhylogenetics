library(dplyr)
library(gghalves)
library(ggnewscale)
library(ggplot2)
library(ggtree)
library(plotly)
library(RColorBrewer)

output <- "/Users/levir/Documents/GitHub/PCAPhylogenetics/manuscript/figures/"

treeSubset <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/sampledTrees.tsv")

#convert to phylo object
phyloList <- list()
for(i in 1:nrow(treeSubset)){
  phyloList[[i]] <- ape::read.tree(text = as.character(treeSubset[i, 2]))
}

# Figure: LDDMM continuous traits --------------------------------------------

tree <- phyloList[[1]]
tree$tip.label <- gsub("_", tree$tip.label, replacement = " ")
p1 <- ggtree(tree)+
  geom_tiplab(fontface = 4)+
  xlim(NA, 5)
p1 <- ggtree::rotate(p1, 36)
p1
ggsave(paste(output, "tree1LDDMMContTraits.svg", sep = ""), p1)

#LDDMMM stacks
lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.100000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
stackFunc <- function(tip) {
  dat <- filter(lddmmRes1Alpha0.2, taxon == tip)
  
  #center dat
  dat$x <- dat$x - mean(dat$x)
  dat$y <- dat$y - mean(dat$y)
  
  y_pad <- 0.5 # Add vertical padding
  y_min <- min(dat$y) - y_pad
  y_max <- max(dat$y) + y_pad
  
  x_pad <- 0.5  # Add vertical padding
  x_min <- min(dat$x) - x_pad
  x_max <- max(dat$x) + x_pad
  
  connections_sim <- dat %>%
    arrange(id) %>%
    mutate(lm_next = lead(id),
           x_next = lead(x),
           y_next = lead(y))
  
  # For the last point, set next to the first point (to close the loop)
  connections_sim[nrow(connections_sim), c("id_next", "x_next", "y_next")] <- 
    c(dat$id[1], dat$x[1], dat$y[1])
  
  ggplot(dat, aes(x = x, y = y)) +
    geom_point(size = 1.5, alpha = 0.8) +
    geom_segment(data = connections_sim,
                 aes(x = x, y = y, xend = x_next, yend = y_next),
                 color = "black", size = 1)+
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme_minimal() +
    ggtitle(tip) +
    theme(
      plot.title = element_text(hjust = 0, size = 9),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
}
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.1.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)


lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.200000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.2.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)

lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.300000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.3.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)

lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.400000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

# Generate plots with fixed axis limits
x_range <- range(lddmmRes1Alpha0.2$x)

# Make individual plots with custom y range
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(paste(output, "fig1LDDMMAlpha0.4.svg", sep = ""), 
       stack_plot, 
       width = 2,
       height = 42)


### 3D
lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/threeDimensionSimTreeIndex0LM10Alpha0.100000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y", "z")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)


# Make individual plots with custom y range
compute_global_xrange <- function(data, taxa) {
  all_proj_x <- c()
  
  for (tip in taxa) {
    dat <- dplyr::filter(data, taxon == tip)
    dat$x <- dat$x - mean(dat$x)
    dat$y <- dat$y - mean(dat$y)
    dat$z <- dat$z - mean(dat$z)
    
    x2d <- dat$x - dat$y * 0.5  # same as in your project_points()
    all_proj_x <- c(all_proj_x, x2d)
  }
  
  return(range(all_proj_x))
}
stackFunc <- function(tip) {
  dat <- filter(lddmmRes1Alpha0.2, taxon == tip)
  
  #center dat
  dat$x <- dat$x - mean(dat$x)
  dat$y <- dat$y - mean(dat$y)
  dat$z <- dat$z - mean(dat$z)
  
  norms <- sqrt(rowSums(dat[,3:5]^2))
  # Divide each row by its norm
  dat[,3:5] <- dat[,3:5] / norms
  
  
  y_pad <- 1.0  # Add vertical padding
  y_min <- min(dat$y) - y_pad
  y_max <- max(dat$y) + y_pad
  
  x_pad <- 1.0  # Add vertical padding
  x_min <- min(dat$x) - x_pad
  x_max <- max(dat$x) + x_pad
  
  
  df <- dat[,3:5]
  
  hull_faces <- geometry::convhulln(df, output.options = TRUE)$hull
  
  get_edges <- function(faces) {
    edges <- do.call(rbind, lapply(1:nrow(faces), function(i) {
      tri <- faces[i, ]
      rbind(tri[c(1, 2)], tri[c(2, 3)], tri[c(3, 1)])
    }))
    edges <- unique(t(apply(edges, 1, sort)))
    return(edges)
  }
  edges <- get_edges(hull_faces)
  hull_vertices <- unique(as.vector(hull_faces))
  
  project_points <- function(x, y, z) {
    x2d <- x - y * 0.5
    y2d <- z - y * 0.5
    return(data.frame(x2d = x2d, y2d = y2d))
  }
  
  projected_df <- project_points(df$x, df$y, df$z)
  
  edge_df <- do.call(rbind, lapply(1:nrow(edges), function(i) {
    idx1 <- edges[i, 1]
    idx2 <- edges[i, 2]
    
    p1 <- df[idx1, ]
    p2 <- df[idx2, ]
    
    proj1 <- projected_df[idx1, ]
    proj2 <- projected_df[idx2, ]
    
    mean_z <- mean(c(p1$z, p2$z))
    
    data.frame(
      x = proj1$x2d, y = proj1$y2d,
      xend = proj2$x2d, yend = proj2$y2d,
      depth = mean_z
    )
  }))
  
  depth_range <- range(edge_df$depth)
  edge_df$gray <- (edge_df$depth - depth_range[1]) / diff(depth_range)
  edge_df$color <- gray(1 - edge_df$gray)
  
  hull_points <- df[hull_vertices, ]
  projected_points <- project_points(hull_points$x, hull_points$y, hull_points$z)
  projected_points$z <- hull_points$z
  projected_points$gray <- (projected_points$z - min(df$z)) / diff(range(df$z))
  projected_points$color <- gray(1 - projected_points$gray)
  
  projected_points <- projected_points %>%
    arrange(z)
  ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x, 
                     y = y, 
                     xend = xend, 
                     yend = yend, 
                     color = depth,
                     alpha = depth),
                 linewidth = 1) +
    scale_alpha_continuous(range = c(0.4, 0.6))+
    geom_point(data = projected_points,
               aes(x = x2d, 
                   y = y2d, 
                   color = z, 
                   size = z)) +
    scale_color_gradient(low = "grey", high = "black") +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme_minimal() +
    ggtitle(tip) +
    theme(
      plot.title = element_text(hjust = 0, size = 9),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position="none",
      plot.margin = margin(2, 2, 2, 2)
    )
}

taxa <- unique(lddmmRes1Alpha0.2$taxon)
x_range <- compute_global_xrange(lddmmRes1Alpha0.2, taxa)
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(filename = paste(output, "3Dfig1LDDMMAlpha01.svg", sep = ""), 
      plot = stack_plot, 
       width = 2,
       height = 42)



lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/threeDimensionSimTreeIndex0LM10Alpha0.200000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y", "z")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

taxa <- unique(lddmmRes1Alpha0.2$taxon)
x_range <- compute_global_xrange(lddmmRes1Alpha0.2, taxa)
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(filename = paste(output, "3Dfig1LDDMMAlpha02.svg", sep = ""), 
       plot = stack_plot, 
       width = 2,
       height = 42)


lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/threeDimensionSimTreeIndex0LM10Alpha0.300000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y", "z")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

taxa <- unique(lddmmRes1Alpha0.2$taxon)
x_range <- compute_global_xrange(lddmmRes1Alpha0.2, taxa)
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(filename = paste(output, "3Dfig1LDDMMAlpha03.svg", sep = ""), 
       plot = stack_plot, 
       width = 2,
       height = 42)


lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/threeDimensionSimTreeIndex0LM10Alpha0.400000Dataset1nodeShapes.tsv", header = FALSE)
colnames(lddmmRes1Alpha0.2) <- c("taxon", "id", "x", "y", "z")
lddmmRes1Alpha0.2$taxon <- gsub("_", lddmmRes1Alpha0.2$taxon, replacement = " ")

tip_order <- p1$data %>%
  filter(isTip) %>%
  arrange(-y) %>%
  pull(label)

taxa <- unique(lddmmRes1Alpha0.2$taxon)
x_range <- compute_global_xrange(lddmmRes1Alpha0.2, taxa)
plots <- lapply(tip_order, stackFunc)

# Stack vertically
stack_plot <-patchwork::wrap_plots(plots, ncol = 1)
stack_plot
ggsave(filename = paste(output, "3Dfig1LDDMMAlpha04.svg", sep = ""), 
       plot = stack_plot, 
       width = 2,
       height = 42)

# Figure: BM continuous traits --------------------------------------------

simResPC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults.rds")
simRes2550PC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults25char50char.rds")
simRes500PC1PC2 <- readRDS("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/simulationResults500char.rds")

simResPC1PC2 <- c(simResPC1PC2, simRes2550PC1PC2, simRes500PC1PC2)

for(i in 1:length(simResPC1PC2)){
  if(i == 1){
    concatDF <- simResPC1PC2[[1]]$resMat
  }else{
    concatDF <- rbind(concatDF, simResPC1PC2[[i]]$resMat)
  }
}
concatDF$numCharacters <- as.factor(concatDF$numCharacters)
concatDF$setRate <- as.factor(concatDF$setRate)
concatDF$propConflicting <- as.factor(concatDF$propConflicting)
concatDF$varRateExpectation <- concatDF$variableRateShape / concatDF$variableRateScale
concatDF$varRateExpectation <- as.factor(concatDF$varRateExpectation)

## set rate vs. variable rates
n_colors <- length(unique(concatDF$setRate))  # Adjust based on your data
blue_colors <- brewer.pal(max(3, n_colors + 2), "Blues")[(3):(n_colors + 2)]  # Skip first 2 colors
red_colors <- brewer.pal(max(3, n_colors + 2), "Reds")[(3):(n_colors + 2)]  # Skip first 2 colors

p1 <- ggplot() +
  geom_half_violin(data = filter(concatDF, is.na(varRateExpectation)),
                  aes(x = numCharacters, y = SPR, fill = setRate), side = "r", scale = "width")+
  scale_fill_manual(name = "Set rate", values = blue_colors) +
  new_scale_fill()+
  geom_half_violin(data = filter(concatDF, is.na(setRate)),
                  aes(x = numCharacters, y = SPR, fill = varRateExpectation), side = "l", scale = "width")+
  scale_fill_manual(name = "Variable rates", values = red_colors) +
  new_scale_fill()+
  scale_y_continuous(breaks = 0:12) +
  xlab(NULL)+
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())
p1

n_colors <- length(unique(concatDF$propConflicting))
blue_colors <- brewer.pal(max(3, n_colors + 2), "Blues")[(3):(n_colors + 2)] 
p2 <- ggplot() +
  geom_violin(data = concatDF,
                  aes(x = numCharacters, y = SPR, fill = propConflicting), scale ="width")+
  scale_fill_manual(name = "% traits conflicting", values = blue_colors) +
  scale_y_continuous(breaks = 0:12) +
  xlab("Number of characters")+
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())
p2

p3 <- p1 / p2
p3

ggsave(paste(output, "varRateSetRateConflictingFigurePC1PC2.svg", sep = ""), p3, width = 10, height = 10)
