# Install & import dependencies ----

#install.packages("tidyverse")
# To get bindrows, %>% and select 
library(dplyr)
library("ape")

install.packages("devtools")
library(devtools)

install_github("vqv/ggbiplot")
#library(ggbiplot)

# Read all the CSV files ----

setwd("~/src/Masters/Project/Notes/Notebooks/Objective 2/data")

file_list <- list.files(pattern="*.pack.table")

read_and_inc_filename <- function(filename) {
  f <- read.delim(filename)
  f$sample <- unlist(strsplit(filename, "[.]"))[[1]]
  f
}

myfiles <- lapply(file_list, read_and_inc_filename)

# Aggregate and normalize ----
aggregate_and_normalize <- function(e) {
  cov <- e[,c("node.id","coverage")]
  
  # Aggregate the median coverage per node
  cov_median <- setNames(aggregate(cov$coverage, list(cov$node.id), median),
                         c("Node ID", "Median Coverage"))
  
  # Extract the sample column to re-append to the output
  sample <- e$sample[1:length(cov_median)]
  
  # Normalize the data 
  total_cov <- sum(cov_median$`Median Coverage`)
  
  normalized_median_coverage <- apply(cov_median[2],
                                      2,
                                      function(x) x/total_cov)
  
  # Return a data frame of the node id against normalized coverage 
  setNames(data.frame(cov_median$`Node ID`, normalized_median_coverage, sample),
                c("Node ID", "Normalized Coverage", "Sample"))
 
}

normalized <- lapply(myfiles, aggregate_and_normalize)

# Binarize normalization
binary_normalization <- function(e) {
  binary_coverage <- apply(e, 1, function(x) if (x==0) {0} else {1})
  
  sample.name <- e$sample[1]
  
  # Return a data frame of the node id against normalized coverage 
  setNames(data.frame(binary_coverage),
           c(sample.name))
}

cov.binary <- lapply(myfiles, binary_normalization)
cov.binary.df <- t(data.frame(cov.binary))

# Tests ----
## Normalization works
assertthat::are_equal(
  length(normalized), 
  sum(unlist(lapply(normalized, function(n) sum(n$`Normalized Coverage`)))))

# Visualize via Heatmap ----
df <- bind_rows(normalized)
ggplot(df, aes(`Node ID`, `Sample`, fill=`Normalized Coverage`)) +
  geom_tile() +
  scale_fill_gradient(low="black", high="white") +
  labs(title = "Normalized and aggregated Coverage vs Node ID per sample",
       x = "Node ID",
       y = "Sample")


# Reshape the data for PCA. Isolate only the coverage. ----
# Extract Normalized Coverage and associate it with the sample name

extract_coverage <- function (e) {
  setNames(data.frame(e %>% select("Normalized Coverage")),
           c(e$Sample[[1]]))
}
samples <- row.names(coverage)

# extract coverage
coverage.list <- lapply(normalized, extract_coverage)
coverage <- data.frame(coverage.list)
coverage <- t(coverage)

# Each household to have its own shape
household <- row.names(coverage)
samples <- row.names(coverage)
# PCA ----
coverage.pca <- prcomp(coverage)
coverage.pca.df <- as.data.frame(coverage.pca$x)

colfunc <- colorRampPalette(c("red", "yellow"))
phage.colors=c(colfunc(10), rainbow(8)[3:43])

library("ggplot2")

# Plot PCA1 against PCA2
ggplot(data = coverage.pca.df, aes(x=PC1, y=PC2, colour=samples)) +
  geom_point() +
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title = "PC1 vs PC2 for Household 20 data")

# Dendrogram/phylogenetics ----
## Install the phylogeny ape package
# install.packages("ape")
install.packages("phyclust")
library(phyclust)

d <- dist(coverage.pca.df)
tree <- nj(d)
plot(tree)
plotnj(tree, show.tip.label = TRUE)
