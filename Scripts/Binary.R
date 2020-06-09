# Imports ----
library("tidyverse")
library("ggplot2")
library(ape)
library(phyclust)


# Utils ----
extract_sample_name <- function(filename) {
  j <- strsplit(filename, "[.]")
  unlist(j)[[1]] 
}


# Extract coverage from the list element and store it in a dataframe
extract_coverage <- function(d) {
  setNames(data.frame(d$coverage),c("coverage")) 
}

# Binarize coverage
binarize_coverage <- function(e) {
  cov.bin <- select(e, .data$coverage) %>% transmute(cov.bin = if_else(coverage==0, 0, 1))
  cov.bin
}

sum_norm_coverage <- function(e) {
  total_cov <- sum(e$coverage)
  select(e, .data$coverage) %>% transmute(cov.norm = coverage/total_cov)
}

sum_norm_bin_coverage <- function(e) {
  total_cov <- sum(e$coverage)
  select(e, .data$coverage) %>% transmute(
    cov.norm.bin = if_else(coverage==0, 0, coverage/total_cov) 
  )
}

# Take f.bin and add a corresponding col to each col with the name of the source sample
append_sample_name_to_bin_cov <- function(name) {
  cov <- select(f.bin, !!name)
  size <- dim(cov)[1]
  
  names <- replicate(size, name)
  positions <- 1:size
  
  setNames(data.frame(positions, cov, names),
           c("Positions","Coverage", "Sample"))
}

# Take f.bin and add a corresponding col to each col with the name of the source sample
append_sample_name_to_bin_cov <- function(name) {
  cov <- select(f.bin, !!name)
  size <- dim(cov)[1]
  
  names <- replicate(size, name)
  positions <- 1:size
  
  setNames(data.frame(positions, cov, names),
           c("Positions","Coverage", "Sample"))
}

append_sample_name_to_named_cov <-  function(name) {
  cov <- x[[name]]
  size <- dim(cov)[1]
  
  names <- replicate(size, name)
  positions <- 1:size
  
  setNames(data.frame(positions, cov, names),
           c("Positions","Coverage", "Sample"))
}

# Script ---- 
## Fetch data ----
setwd("~/src/Masters/Project/Notes/Notebooks/Objective 2/data")

# Get the list of each coverage table
file_list <- list.files(pattern="*.pack.table")

# Extract each sample name
sample_names <- lapply(file_list, extract_sample_name)

# read all the CSV files into a list
all_columns <- lapply(file_list, read.delim)

# Extract the coverage column for each sample
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
ggsave("~/src/Masters/Project/Notes/Images/Objective 2/RSV/Median_Coverage_HH20.png", 
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
ggsave("~/src/Masters/Project/Notes/Images/Objective 2/RSV/Mean_Coverage_HH20.png", 
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
ggsave("~/src/Masters/Project/Notes/Images/Objective 2/RSV/Max_Coverage_HH20.png", 
       height=6, 
       width=8, 
       dpi=400)
dev.off()

ggplot(data=covs.structure) +
  geom_col(aes(x=Samples, y=variance)) +
  theme(axis.text.x = element_text(angle = 90))+
  labs(title = "Standard Deviation of Coverage per sample of the dataset",
       x = "Sample",
       y = "Standard Deviation of Sample Coverage")
ggsave("~/src/Masters/Project/Notes/Images/Objective 2/RSV/Std_Coverage_HH20.png", 
       height=6, 
       width=8, 
       dpi=400)
dev.off()

## Normalization ----

## Binary normalize the coverage
### if coverage value is above 0 convert it to 1
f.bin <- setNames(data.frame(lapply(covs, binarize_coverage)), sample_names)
ft.bin <- t(f.bin)

## Normalize so the the total coverage of each sample is 1
f.norm <- setNames(data.frame(lapply(covs, sum_norm_coverage)), sample_names)
ft.norm <- t(f.norm)

## sum normalize then binarize
f.norm.bin <- setNames(data.frame(lapply(covs, sum_norm_bin_coverage)), sample_names)
ft.norm.bin <- t(f.norm.bin)

## Heatmap ----
### SARS-CoV-2 ---
x <- setNames(covs, sample_names)
x.named <- lapply(sample_names, append_sample_name_to_named_cov)
bound_samples <- bind_rows(x.named)

### RSV ---
lp <- lapply(sample_names, append_sample_name_to_bin_cov)
bound_samples <- bind_rows(lp)

## Visualize ---
ggplot(bound_samples, aes(`Positions`, `Sample`, fill=`Coverage`)) +
  geom_tile() +
  scale_fill_gradient(low="black", high="white") +
  labs(title = "SARS-COV-2 simulated reads Coverage vs Position per sample",
       x = "Position",
       y = "Sample")
ggsave("~/src/Masters/Project/Notes/Images/Objective 2/SARS_Heatmap.png", 
       height=6, 
       width=8, 
       dpi=400)
dev.off()

## PCA ----
### SARS-CoV-2

xt.df.pca <- prcomp(t(setNames(data.frame(x), sample_names)))
data <- as.data.frame(xt.df.pca$x) 
  
### RSV
ft.bin.pca <- prcomp(ft.bin)
ft.bin.pca.df <- as.data.frame(ft.bin.pca$x)
data <- ft.bin.pca.df

colfunc <- colorRampPalette(c("red", "yellow"))
phage.colors=c(colfunc(10), rainbow(8)[3:43])

# Plot PCA1 against PCA2
Samples <- unlist(sample_names)

ggplot(data = data, aes(x=PC1, y=PC2, colour=Samples)) +
  geom_point() +
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title = "SARS-CoV-2 simulated reads PC1 vs PC2 ")
ggsave("~/src/Masters/Project/Notes/Images/Objective 2/SARS_CoV_2_PCA.pdf", 
       height=6, 
       width=8, 
       dpi=400)
dev.off()

## Neighbour joining ----

unique_households <- unique(lapply(sample_names, function(n) { unlist(strsplit(n[[1]], "_"))[[2]]}))

d <- dist(ft.bin.pca.df)
tree <- nj(d)

write.tree(tree, file="~/Desktop/my-tree.nwk", tree.names = TRUE)

pdf(file="~/src/Masters/Project/Notes/Images/Objective 2/household_20_rooted_dendrogram.pdf", 
    height=18, 
    width=18)
colors <- sample(colors(), 10)
plot(tree, 
     "u",
     # cex=2, # ?tiplabels
     # adj = c(0, 0.5),
     main = "Phylogenetic tree for household 20")

dev.off()
