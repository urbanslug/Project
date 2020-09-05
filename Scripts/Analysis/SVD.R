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
output_dir <- "~/src/Masters/Project/Figures/SVD/"
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
individuals <- unlist(lapply(household_members, extract_individuals))

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
colnames(coverage.df.orig) <- household_members

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




# TODO: what is it really doing with these many cols?
coverage.rsvd.large <- rsvd(coverage.matrix, k=cols)
lv <- coverage.rsvd.large$v
lv.df <- data.frame(lv, row.names = metadata.gg$samples)
colnames(lv.df) <- metadata.gg$samples
# lv.df$household <- factor(metadata.gg$households)

# Visualiztion ----
## 
lv.m <- as.dist(lv.df)
lv.m.new <- dist(lv.m)
lv.tree.m <- nj(lv.m.new)

write.tree(lv.tree.m, "~/src/Masters/Project/Trees/rsv_svd_tree.nwk")

myPalette <-c("#708090", "#0014a8", "#9f00ff", "#177245", "#f984ef", "#ffae42", 
              "#03c03c", "#915f6d", "#f7e98e", "#0070ff", "#663854", "#e8000d",
              "#704214", "#00ced1", "#ffa07a", "#b5651d",  "#918151")

p <- ggtree(lv.tree.m) %<+% metadata.gg
p +
  geom_tiplab(size=4) + 
  geom_tippoint(size=2, aes(color=individuals)) +
  labs(title = "RSV Neighbour Joining Tree of Individuals Sampled Over Time", 
       color="Individuals") +
  scale_color_manual(values=myPalette) +
  theme_tree(legend.position='left')

ggsave(paste(output_dir, "RSV_SVD_nj_Tree.png", sep=""),
       height=10, 
       width=22, 
       dpi=300)
dev.off()


# PCA
coverage.pca<- rpca(t(coverage.matrix), k=cols)
coverage.pca.x <- coverage.pca$x
coverage.pca.dist <- dist(coverage.pca.x)
coverage.pca.tree <- nj(coverage.pca.dist)

write.tree(coverage.pca.tree, "~/src/Masters/Project/Trees/RSV_PCA_Tree.nwk")
write.csv(metadata.gg, "~/src/Masters/Project/Data/RSV_Metadata.csv")

p <- ggtree(coverage.pca.tree) %<+% metadata.gg
p +
  geom_tiplab(size=4) + 
  geom_tippoint(size=2, aes(color=individuals)) +
  labs(title = "RSV PCA Neighbour Joining Tree", color="Individuals") +
  scale_color_manual(values=myPalette) +
  theme_tree(legend.position='left')

ggsave(paste(output_dir, "RSV_PCA_nj_Tree_29.png", sep=""),
       height=10, 
       width=22, 
       dpi=300)
dev.off()


# Hierarchical clustering -----
x <- hclust(lv.m.new)
plot(x)

# Please ignore ----------------------------------------------------------------
# Vary number of SVD 
coverage.rsvd <- rsvd(coverage.matrix, k=29)
sv <- coverage.rsvd$v
sv.df <- data.frame(sv, row.names = metadata.gg$samples)

# sv.dist_matrix <- as.dist(sv.df)
sv.m <- as.dist(sv.df)
sv.m.new <- dist(sv.m)
sv.tree.m <- nj(sv.m.new)

p <- ggtree(sv.tree.m) %<+% metadata.gg
p +
  geom_tiplab(size=4) + 
  geom_tippoint(size=2, aes(color=individuals)) +
  labs(title = "RSV SVD Neighbour Joining Tree", color="Individuals") +
  scale_color_manual(values=myPalette) +
  theme_tree(legend.position='left')