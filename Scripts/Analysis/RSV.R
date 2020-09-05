#!/usr/bin/env Rscript

# Imports ----
suppressPackageStartupMessages(
  {
    require(dplyr)
    require(FactoMineR)
    require(ggplot2)
    require(ggtree)
    require(ape)
  })

# Global Vars ----
output_dir <- "~/src/Masters/Project/Figures/RSV"

## Fetch data ----
setwd("~/src/Masters/Project/Notebooks/Objective 2/data")

# Get the list of filenames ending in .pack.table
file_list <- list.files(pattern="*.pack.table")

# read all the CSV files into a list
all_columns <- lapply(file_list, read.delim)

# Extract the coverage column for each sample
# Extract coverage from the list element and store it in a dataframe
extract_coverage <- function(d) {
  setNames(data.frame(d$coverage),c("coverage")) 
}
covs <- lapply(all_columns, extract_coverage)



## Structure of the data ----
medians <- as.vector(unlist(lapply(covs, function(x) {median(x[,])})))
means <- as.vector(unlist(lapply(covs, function(x) {mean(x[,])}))) 
maxs <- as.vector(unlist(lapply(covs, function(x) {max(x[,])}))) 

sd <- as.vector(unlist(lapply(covs, function(x) {sd(x[,])}))) 

covs.structure <- data.frame(medians, means, maxs, sd, c(unlist(sample_names)))
colnames(covs.structure) <- c("Median", "Mean", "Max", "Std Dev",  "Samples")
row.names(covs.structure) <- sample_names

### Visualize the structure of the data
ggplot(data=covs.structure) +
  geom_col(aes(x=Samples, y=medians)) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title = "Meadian Coverage per sample of the dataset",
       x = "Sample",
       y = "Median Sample Coverage")
ggsave(paste(output_dir, "Median_Coverage_HH20.png", sep=""), 
       height=6, 
       width=8, 
       dpi=400)
dev.off()

ggplot(data=covs.structure) +
  geom_col(aes(x=Samples, y=means)) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title = "Mean Coverage per sample of the dataset",
       x = "Sample",
       y = "Mean Sample Coverage")
ggsave("~/src/Masters/Project/Figures/RSV/Mean_Coverage_HH20.png", 
       height=6, 
       width=8, 
       dpi=400)
dev.off()

ggplot(data=covs.structure) +
  geom_col(aes(x=Samples, y=maxs)) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title = "Max Coverage per sample of the dataset",
       x = "Sample",
       y = "Maximum coverage value in the sample")
ggsave("~/src/Masters/Project/Figures/RSV/Max_Coverage_HH20.png", 
       height=6, 
       width=8, 
       dpi=400)
dev.off()

ggplot(data=covs.structure) +
  geom_col(aes(x=Samples, y=sd)) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title = "Standard Deviation of Coverage per sample of the dataset",
       x = "Sample",
       y = "Standard Deviation of Sample Coverage")
ggsave("~/src/Masters/Project/Figures/RSV/Std_Coverage_HH20.png", 
       height=6, 
       width=8, 
       dpi=400)
dev.off()


## Normalization ----

## Binary normalize the coverage
### if coverage value is above 0 convert it to 1

# Extract each sample name
extract_sample_name <- function(filename) {
  j <- strsplit(filename, "[.]")
  unlist(j)[[1]] 
}
sample_names <- lapply(file_list, extract_sample_name)

# Binarize coverage
binarize_coverage <- function(e) {
  cov.bin <- select(e, .data$coverage) %>% transmute(cov.bin = if_else(coverage==0, 0, 1))
  cov.bin
}

f.bin <- setNames(data.frame(lapply(covs, binarize_coverage)), sample_names)
ft.bin <- t(f.bin)

# PCA ----
ft.bin.pca <- PCA(ft.bin)

# Phylogenetics ----
## Extract the dimensions of the PCA
individuals <- ft.bin.pca$ind
dimensions <- data.frame(individuals[["coord"]])

## Prepare metadata ----
households <- lapply(sample_names, function(n) { unlist(strsplit(n[[1]], "_"))[[2]]})
households.df <- setNames(do.call(rbind.data.frame, households), c("hh"))
l <- setNames(data.frame(cbind(unlist(sample_names), unlist(households))), 
         c("sample_names", "Households"))

## First principal component ----
pc1 <- dimensions[1]
m.pc1 <- as.matrix(pc1)
distance.matrix.pc1 <- dist(m.pc1)
tree.pc1 <- nj(distance.matrix.pc1)

p <- ggtree(tree.pc1) %<+% l
p + geom_treescale() +
  geom_tiplab(size=4, aes(color=Households)) + 
  geom_tippoint(size=1, aes(color=Households)) +
  labs(title = "Neighbour Joining Tree Based on Graph Coverage and First Principal Component") +
  theme_tree(legend.position=c(.9, .3))

ggsave("~/src/Masters/Project/Figures/RSV/HH20_first_principal_component_nj_tree.png", 
       height=10, 
       width=22, 
       dpi=300)
dev.off()

## All principal components ----
m.all_pcs <- as.matrix(dimensions)
distance.matrix.all_pcs <- dist(m.all_pcs)
tree.all_pcs <- nj(distance.matrix)

write.tree(tree.all_pcs, file="~/Desktop/all_pcs.nwk", tree.names = TRUE)

a <- ggtree(tree.all_pcs) %<+% l
a + geom_treescale() +
  geom_tiplab(size=4, aes(color=Households)) + 
  geom_tippoint(size=1, aes(color=Households)) +
  labs(title = "Neighbour Joining Tree Based on Graph Coverage and All Principal Components") +
  theme_tree(legend.position=c(.9, .3))

ggsave("~/src/Masters/Project/Figures/RSV/HH20_all_principal_components_nj_tree.png", 
       height=10, 
       width=22, 
       dpi=300)
dev.off()

## Old tree ----

old_tree <- read.tree("~/Desktop/Trees/my-tree.nwk")
ggtree(old_tree) %<+% l +
  geom_treescale() + 
  geom_tiplab(aes(color=Households)) +
  geom_tippoint(aes(color=Households)) +
  labs(title = "RSV ")
