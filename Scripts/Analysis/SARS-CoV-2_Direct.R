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
### SARS-CoV-2
output_dir <- "~/src/Masters/Project/Figures/Analysis/SARS-CoV-2/"
setwd("~/src/Masters/Project/Data/SARS-CoV-2/")

# Fetch data ----
# Get the list of filenames ending in .pack.table
file_list <- list.files(pattern="*.pack.table")

# read all the CSV files into a list
all_columns <- lapply(file_list, read.delim)

# Extract the coverage column for each sample
# Extract coverage from the list element and store it in a dataframe
extract_coverage <- function(d) {
  setNames(data.frame(d$coverage),c("coverage")) 
}
coverage.l <- lapply(all_columns, extract_coverage)

# Mangle it up ----
# Extract each sample name 
extract_sample_name <- function(filename) {
  j <- strsplit(filename, "[.]")
  n <- unlist(j)[[1]] 
  sample_num <- unlist(strsplit(n, "-"))[4]
  paste("Sample", sample_num, sep=" ")
}
sample_names <- unlist(lapply(file_list, extract_sample_name))
samples_names.df <- data.frame(sample_names)

binarize_coverage <- function(e) {
  # Given a list of coverage values convert it to a vector of zeros and ones
  binarized <- e %>% transmute(cov.bin = if_else(coverage==0, 0, 1))
  binarized
}

coverage.binary <- lapply(coverage.l, binarize_coverage)

coverage.df <- data.frame(coverage.binary)
colnames(coverage.df) <- sample_names

coverage.matrix <- as.matrix(coverage.df)
rows <- dim(coverage.matrix)[1]
cols <- dim(coverage.matrix)[2]

# Direct ----
coverage.matrix.dist <-dist(t(coverage.matrix))

## Neighbor Joining ----
coverage.nj <- nj(coverage.matrix.dist)
coverage.nj.ggtree <- ggtree(coverage.nj)  %<+% samples_names.df
coverage.nj.ggtree + 
  geom_tiplab(size=3, offset=0) + 
  geom_tippoint(size=1, aes(color=label)) +
  labs(title = "SARS-CoV-2 Neighbour Joining Cladogram", color="Samples") +
  theme_tree(legend.position='')

ggsave(paste(output_dir, "SARS-CoV-2_NJ.png", sep=""),
       height=5, 
       width=12, 
       dpi=400)
dev.off()

## Hierachical Clustering ----
coverage.clust <- hclust(coverage.matrix.dist)
coverage.hclust.ggtree <- ggtree(coverage.clust) %<+% samples_names.df
coverage.hclust.ggtree + geom_tiplab(size=5, offset=0) + 
  geom_tippoint(size=3, aes(color=label)) +
  labs(title = "SARS-CoV-2 Hierachical Clustering", color="Samples") +
  theme_tree(legend.position='')

ggsave(paste(output_dir, "SARS-CoV-2_HClust.png", sep=""),
       height=8, 
       width=20,
       units="in",
       limitsize = FALSE,
       dpi=200)
dev.off()