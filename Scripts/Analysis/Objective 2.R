# Install dependencies
install.packages("tidyverse")

# Read all the CSV files
setwd("~/src/Masters/Project/Notes/Notebooks/Objective 2/data")

file_list <- list.files(pattern="*.pack.table")

read_and_inc_filename <- function(filename) {
  f <- read.delim(filename)
  f$sample <- unlist(strsplit(filename, "[.]"))[[1]]
  f
}

myfiles <- lapply(file_list, read_and_inc_filename)



aggregate_and_normalize <- function(e) {
  cov <- e[,c("node.id","coverage")]
  
  # Aggregate the median coverage per node
  cov_median <- setNames(aggregate(cov$coverage, list(cov$node.id), median),
                         c("Node ID", "Median Coverage"))
  
  
  # Normalize the data 
  total_cov <- sum(cov_median$`Median Coverage`)
  
  normalized_median_coverage <- apply(cov_median[2],
                                      2,
                                      function(x) x/total_cov)
  
  # Return a data frame of the node id against normalized coverage 
  f <- setNames(data.frame(cov_median$`Node ID`, normalized_median_coverage),
           c("Node ID", "Normalized Coverage"))
  f$sample <- e$sample[1:length(cov_median)]
  f
}

normalized <- lapply(myfiles, aggregate_and_normalize)

library(dplyr)

df <- bind_rows(normalized)

ggplot(df, aes(`Node ID`, `sample`, fill=`Normalized Coverage`)) +
  geom_tile() +
  scale_fill_gradient(low="black", high="white") +
  labs(title = "Coverage...",
       x = "Node ID",
       y = "Sample")

norm_len <- length(normalized)

normalized 

function

samples <- lapply(file_list, function(name) replicate(length(normalized), name))



lapply(normalized, function(x) x)


df <- do.call(rbind, lapply(file_list, function(x) read.csv(x, sep="\t")))


library(tidyverse)



# Associate each ... with the file in the file list
attach_names <- function(s) {
  l <- dim(s)
  replicate(node.num, s)
  
} 

tbl <-
  list.files(pattern = "*.pack.table") %>% 
  map_df(~read_csv(.))


setwd("~/src/Masters/Project/Notes/Scripts")
H_503_16_04 <- read.csv("../Notebooks/Objective 2/data/H_503_16_04.pack.table", 
                        sep="\t")
H_502_09_04 <- read.csv("../Notebooks/Objective 2/data/H_502_09_04.pack.table", 
                        c)

cov <- data.frame(H_503_16_04[,c("node.id","coverage")])
cov2 <- data.frame(H_502_09_04[,c("node.id","coverage")])

# Get the median coverage for each node ID
median_coverage <- setNames(aggregate(cov$coverage, list(node_id=cov$node.id), median), 
                            c("Node ID", "Median Coverage"))
median_coverage2 <- setNames(aggregate(cov2$coverage, list(node_id=cov2$node.id), median), 
                            c("Node ID", "Median Coverage"))

typeof(median_coverage)


# View(median_coverage)
med_cov <- median_coverage$`Median Coverage`

# View(med_cov)

# Type of med_cov is a list of median coverage
med_cov_max <- max(median_coverage$`Median Coverage`)
med_cov_min <- min(median_coverage$`Median Coverage`)
med_cov_sum <- sum(median_coverage$`Median Coverage`)
med_cov_sum2 <- sum(median_coverage2$`Median Coverage`)


normalized_cov <- apply(median_coverage[2],
                        2,
                        function(x) x/med_cov_sum)
normalized_cov2 <- apply(median_coverage2[2],
                        2,
                        function(x) x/med_cov_sum2)

node.num <- length(normalized_cov)

l <- 1:node.num

a <- setNames(data.frame(l, normalized_cov, replicate(node.num, "H_503_16_04")), 
              c("Node ID","Coverage","Sample")) 

b <- setNames(data.frame(1:length(normalized_cov2), normalized_cov2, replicate(length(normalized_cov2), "H_502_09_04")),
              c("Node ID","Coverage","Sample")) 

df <- rbind(a, b)

ggplot(df, aes(`Node ID`, Sample, fill=Coverage)) +
  geom_tile()+
  scale_fill_gradient(low="black", high="white") +
  # scale_fill_distiller(palette = "Greys") +
  labs(title = "Coverage",
       x = "Node ID",
       y = "Sample")

sum(normalized_cov)
sum(normalized_cov2)

typeof(normalized_cov)
View(normalized_cov)
View(normalized_cov2)

x<-matrix(data=c(normalized_cov, normalized_cov2))

data.frame(t(matrix(normalized_cov)), 

df <- setNames(data.frame(normalized_cov, normalized_cov2), c("H_503_16_04", "H_502_09_04"))
View(df)
s <-data.frame(transpose(df))
View(s)
typeof(s)
normalized_coverage <- setNames(data.frame(median_coverage$`Node ID`, normalized_cov), c("NodeID", "NormalizedCoverage"))
View(normalized_coverage)
typeof(normalized_coverage)



# Visualizations
# ----
## If not installed
## install.packages("ggplot2")

library("ggplot2")

ggplot(data=s) +
  geom_tile()+ 
  scale_fill_distiller(palette = "YlGnBu") +
  labs(title = "Coverage of Household data",
       x = "Node Identifier",
       y = "Normalized Coverage")

View(mpg)

typeof(mpg)

ggplot(data=mpg) +
  geom_point(mapping=aes(y=manufacturer, x=year)) +
  labs(title = "Likelihood of swinging and missing on a fastball",
       x = "manufacturer",
       y = "year")
