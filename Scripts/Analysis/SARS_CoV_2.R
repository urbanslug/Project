#!/usr/bin/env Rscript

# Imports ----
suppressPackageStartupMessages(
  {
    require(dplyr)
    require(FactoMineR)
    require(ggplot2)
    require(ggtree)
  })

# Global Vars ----
output_dir <- "~/src/Masters/Project/Figures/SARS_CoV_2/"

## Fetch data ----
setwd("~/src/Masters/Project/Notebooks/Objective 1/data")

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

# Extract each sample name
extract_sample_name <- function(filename) {
  j <- strsplit(filename, "[.]")
  n <- unlist(j)[[1]] 
  sample_num <- unlist(strsplit(n, "-"))[4]
  paste("Sample", sample_num, sep=" ")
}
sample_names <- lapply(file_list, extract_sample_name)

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

# PCA ----
x <- setNames(covs, sample_names) # a list
x.df <- setNames(data.frame(x), sample_names)

sars.pca <- PCA(t(x.df), graph = FALSE)

## Extract the dimensions of the PCA
individuals <- sars.pca$ind
dimensions.sars <- data.frame(individuals[["coord"]])

# Append the sample names column
dimensions.sars$Samples <- unlist(sample_names)

ggplot(dimensions.sars, aes(x=Dim.1, y=Dim.2, colour=Samples)) +
  geom_point()

# Phylogenetics ----


## Prepare metadata ----
sample_number <- lapply(sample_names, 
                        function(n) { unlist(strsplit(n[[1]], "-"))[[1]]})
l <- setNames(data.frame(cbind(unlist(sample_names), unlist(sample_number))), 
              c("sample_names", "Samples"))

## First principal component ----
pc1 <- dimensions.sars[1]
m.pc1 <- as.matrix(pc1)
distance.matrix.pc1 <- dist(m.pc1)
tree.pc1 <- nj(distance.matrix.pc1)

p <- ggtree(tree.pc1) %<+% l
p + geom_treescale() +
  geom_tiplab(size=4, aes(color=Samples)) + 
  geom_tippoint(size=1, aes(color=Samples)) +
  labs(title = "SARS Cov 2 Neighbour Joining Tree Based on Graph Coverage and First Principal Component") +
  theme_tree(legend.position=c(.9, .7))

ggsave(paste(output_dir, "SARS_CoV_2_first_principal_component_nj_tree.png", sep=""),
       height=10, 
       width=22, 
       dpi=300)
dev.off()

## All principal components ----
sars.all_pcs <- as.matrix(dimensions.sars)
distance.matrix.all_pcs <- dist(sars.all_pcs)
tree.all_pcs.sars <- nj(distance.matrix.all_pcs)

b <- ggtree(tree.all_pcs.sars) %<+% l
b + geom_treescale() +
  geom_tiplab(size=4, aes(color=Samples)) + 
  geom_tippoint(size=1, aes(color=Samples)) +
  labs(title = "SARS Cov 2 Neighbour Joining Tree Based on Graph Coverage and All Principal Components") +
  theme_tree(legend.position=c(.9, .3))

ggsave(paste(output_dir, "SARS_CoV_2_all_principal_components_nj_tree.png", sep=""),
       height=10, 
       width=22, 
       dpi=300)
dev.off()
