#!/usr/bin/env Rscript

# Imports ----
suppressPackageStartupMessages(
  {
    require(dplyr)
    require(FactoMineR)
    require(ggplot2)
    require(ggtree)
    require(rsvd)
    require(FactoMineR)
    require(pracma)
    require(ape)
  })

# Set up ----
## Vars
### RSV
output_dir <- "~/src/Masters/Project/Figures/Analysis/RSV/"
setwd("~/src/Masters/Project/Data/RSV/")


# Fetch data ----
# Get the list of file names in the current dir ending in .pack.table
file_list <- list.files(pattern="*.pack.table")

## Metadata ----
# Create sample names from the file list
# Extract each sample name
extract_samples <- function(filename) {
  j <- strsplit(filename, "[._]")
  k<-j[[1]][2:4]
  paste(k, collapse="_")
}
samples <- unlist(lapply(file_list, extract_samples))

extract_individuals <- function(n) { unlist(strsplit(n[[1]], "_"))[[1]] }
individuals <- unlist(lapply(samples, extract_individuals))

metadata <- data.frame(samples, individuals)

## Actually fetch the data ----
## Read all the CSV files into a list
## TODO: rename all_columns
all_columns <- lapply(file_list, read.delim)


# Extract the coverage column for each sample
# Extract coverage from the list element and store it in a dataframe
extract_coverage <- function(d) {
  setNames(data.frame(d$coverage),c("coverage"))
}

## Extract the coverage column into a list & associate it with its sample name
coverage.l <- lapply(all_columns, extract_coverage)

# Normalization ----
## Binary normalization
## Binarize coverage
binarize_coverage <- function(e) {
  # Given a list of coverage values convert it to a vector of zeros and ones
  binarized <- e %>% transmute(cov.bin = if_else(coverage==0, 0, 1))
  binarized
}

coverage.binary <- lapply(coverage.l, binarize_coverage)

# Align with data from Githinji 2018 ----
coverage.df.orig <- data.frame(coverage.binary)
colnames(coverage.df.orig) <- samples

# These are the individuals from the paper 
# https://doi.org/10.1101/411512 Githinji 2018
ggcolsorig <- c("525_05_04", "525_02_04", "531_13_04", "506_09_04", "511_13_04",
                "509_02_04", "503_16_04", "503_20_04", "508_02_04", "514_21_05",
                "507_09_04", "504_02_04", "536_09_04", "513_30_03", "515_30_03",
                "513_05_04", "513_02_04", "512_23_03", "518_19_03", "518_16_03",
                "506_30_04", "512_26_03", "504_26_03", "529_02_04", "528_13_04",
                "502_09_04", "507_30_04", "515_02_04", "508_09_04", "507_05_04", 
                "506_13_04", "518_09_03")
ggcolsorig <- sort(ggcolsorig)
# Those that we don't have in common
difference <- c("536_09_04", "507_30_04", "506_13_04")
ggcols <- setdiff(ggcolsorig, difference)

coverage.df <- coverage.df.orig[, ggcols]
metadata.gg <- subset(metadata, samples %in% ggcols)

# Analysis ----
coverage.matrix <- as.matrix(coverage.df)
rows <- dim(coverage.matrix)[1]
cols <- dim(coverage.matrix)[2]


# Visualiztion ----
# From https://sashamaps.net/docs/resources/20-colors/
palette <-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
            '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
            '#008080', '#e6beff', '#9a6324', '#ffd8b1', '#000000',
            '#aaffc3', '#808000')

## PCA ----
coverage.pca <- rpca(t(coverage.matrix), k=cols)
coverage.pca.x <- coverage.pca$x
coverage.pca.df <- data.frame(coverage.pca.x)
coverage.pca.df.gg <- cbind(metadata.gg, coverage.pca.df)

### plot the Principal Components
#### PC1 vs PC2
ggplot(data=coverage.pca.df.gg) +
  geom_point(aes(x=X1, y=X2, color=individuals)) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title = "Principal Components 1 vs 2 For RSV HouseHold 5 Datataset",
       x = "PC1", 
       y = "PC2", 
       color="Individuals")

ggsave(paste(output_dir, "RSV_PCA1vs2.png", sep=""), height=10, width=22, dpi=300)
dev.off()

#### PC3 vs PC4
ggplot(data=coverage.pca.df.gg) +
  geom_point(aes(x=X3, y=X4, color=individuals)) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title = "Principal Components 3 vs 4 For RSV HouseHold 5 Datataset",
       x = "PC3", 
       y = "PC4", 
       color="Individuals")

ggsave(paste(output_dir, "RSV_PCA3vs4.png", sep=""), height=10, width=22, dpi=300)
dev.off()

## Hierarchical clustering -----
coverage.pca.dist <- dist(coverage.pca.x)
coverage.pca.dist.scaled <- as.dist(scale(coverage.pca.dist, center=FALSE))
coverage.pca.hclust <- hclust(coverage.pca.dist)
a <- ggtree(coverage.pca.hclust) %<+% metadata.gg
a + geom_tiplab(size=4, offset=0) + 
  geom_tippoint(size=2, aes(color=individuals)) +
  labs(title = "RSV PCA Hierachical Clustering Tree", color="Individuals") +
  scale_color_manual(values=palette) +
  theme_tree(legend.position='left')

ggsave(paste(output_dir, "RSV_PCA_HC_Tree2.png", sep=""),
       height=10, 
       width=22, 
       dpi=300)
dev.off()


X <-dist(t(coverage.matrix))

i <- hclust(X)
j <- ggtree(i) %<+% metadata.gg
j + geom_tiplab(size=4, offset=0) + 
  geom_tippoint(size=2, aes(color=individuals)) +
  labs(title = "RSV PCA Hierachical Clustering Tree", color="Individuals") +
  scale_color_manual(values=palette) +
  theme_tree(legend.position='left')

ggsave(paste(output_dir, "RSV_PCA_HC_Tree3.png", sep=""),
       height=10, 
       width=22, 
       dpi=300)
dev.off()

X.tree <- nj(X)

p <- ggtree(X.tree)  %<+% metadata.gg
p + geom_tiplab(size=4, offset=0) + 
  geom_tippoint(size=2, aes(color=label)) +
  labs(title = "SARS-CoV-2 PCA Neighbour Joining Tree", color="Samples") +
  theme_tree(legend.position='left')

ggsave(paste(output_dir, "RSV_PCA_HC_Tree4.png", sep=""),
       height=10, 
       width=22, 
       dpi=300)
dev.off()