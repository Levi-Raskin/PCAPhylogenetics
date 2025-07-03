<<<<<<< Updated upstream
<<<<<<< Updated upstream
library(dplyr)
library(ggplot2)
library(ggtree)

output <- "/Users/levir/Documents/GitHub/PCAPhylogenetics/manuscript/figures/"
=======
=======
>>>>>>> Stashed changes
library(ggplot2)
library(ggtree)

output <- "/Users/levir/Documents/GitHub/PCAPhylogenetics/manuscript/figures"
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

treeSubset <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/sampledTrees.tsv")

#convert to phylo object
phyloList <- list()
for(i in 1:nrow(treeSubset)){
  phyloList[[i]] <- ape::read.tree(text = as.character(treeSubset[i, 2]))
}


<<<<<<< Updated upstream
<<<<<<< Updated upstream
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
lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.200000Dataset0nodeShapes.tsv", header = FALSE)
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
  
  y_pad <- 0.2  # Add vertical padding
  y_min <- min(dat$y) - y_pad
  y_max <- max(dat$y) + y_pad
  
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
    coord_cartesian(xlim = x_range, ylim = c(y_min, y_max)) +
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
ggsave(paste(output, "fig1LDDMMAlpha0.2.svg", sep = ""), 
       stack_plot, 
       width = 5,
       height = 30)


lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.400000Dataset0nodeShapes.tsv", header = FALSE)
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
       width = 5,
       height = 30)
lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.600000Dataset0nodeShapes.tsv", header = FALSE)
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
ggsave(paste(output, "fig1LDDMMAlpha0.6.svg", sep = ""), 
       stack_plot, 
       width = 5,
       height = 30)

lddmmRes1Alpha0.2 <- read.delim("/Users/levir/Documents/GitHub/PCAPhylogenetics/results/Mongle_et_al_2023_RB/SimRes/LDDMMSimRes/twoDimensionSimTreeIndex0LM10Alpha0.800000Dataset0nodeShapes.tsv", header = FALSE)
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
ggsave(paste(output, "fig1LDDMMAlpha0.8.svg", sep = ""), 
       stack_plot, 
       width = 5,
       height = 30)

=======
=======
>>>>>>> Stashed changes
# Figure: BM continuous traits --------------------------------------------

tree <- phyloList[[1]]
tree$tip.label <- gsub("_", tree$tip.label, replacement = " ")
<<<<<<< Updated upstream
tree <- ggtree()
>>>>>>> Stashed changes
=======
tree <- ggtree()
>>>>>>> Stashed changes
