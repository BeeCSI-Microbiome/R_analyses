library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# read data
ctx_data <- read.delim(file = '2020_ctx_kraken2/ctx_kraken_genus_data.tsv',
                       header = TRUE,
                       sep = '\t')

# clean data
clean_ctx_data <- select(ctx_data, -taxRank, -taxID, -lineage) %>%
  filter(name %in% c("Gilliamella", 
                  "Frischella", 
                  "Snodgrassella", 
                  "Lactobacillus", 
                  "Bifidobacterium")) %>%
  pivot_longer(!name, names_to = "sample", values_to = "count") %>%
  pivot_wider(names_from = "name", values_from = "count")

# add treatment info
treatments <- c("Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI")
clean_ctx_data$treatment <- treatments

# for ordering on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
clean_ctx_data$treatment <- factor(clean_ctx_data$treatment,
                                   levels = order)

# add replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")
clean_ctx_data$replicate <- replicates

# reorder data frame
clean_ctx_data <- clean_ctx_data[, c(1,8,7,2,3,4,5,6)]

# convert genus part of data frame to matrix for nmds
genus_data <- clean_ctx_data[,4:ncol(clean_ctx_data)]
genus_mat <- as.matrix(genus_data) 

# perform nmds
set.seed(1)
nmds = metaMDS(genus_mat, distance = "bray")

# add back information to matrix
nmds_data = as.data.frame(scores(nmds))
nmds_data$sample = clean_ctx_data$sample
nmds_data$replicate = clean_ctx_data$replicate
nmds_data$treatment = clean_ctx_data$treatment

# plot data
nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = treatment, colour = replicate), size = 5) +
  labs(title = "NMDS Analysis",
       x = "NMDS1", 
       y = "NMDS2", 
       shape = "Treatment", 
       colour = "Replicate")

nmds_plot

# save plot as svg
# svg("NMDS_Plot.svg")
# nmds_plot
# dev.off()