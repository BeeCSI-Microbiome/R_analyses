library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# read data
ctx_data <- read.delim(file = 'ctx_kraken_species_precent_data.tsv',
                       header = TRUE,
                       sep = '\t')

# clean data
clean_ctx_data <- select(ctx_data, -taxRank, -taxID, -Max, -lineage) %>%
  filter(name != "Apis mellifera") %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")

# convert taxa part of data frame to matrix for nmds
taxa_data <- clean_ctx_data[,2:ncol(clean_ctx_data)]
taxa_mat <- as.matrix(taxa_data) 

# perform nmds
set.seed(1)
nmds = metaMDS(taxa_mat, distance = "bray")

# add treatment info
treatments <- c("Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI")

# add replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")

# add back information to matrix
nmds_data = as.data.frame(scores(nmds))
nmds_data$sample = clean_ctx_data$sample
nmds_data$replicate <- replicates
nmds_data$treatment <- treatments

# for ordering on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
nmds_data$treatment <- factor(nmds_data$treatment,
                              levels = order)

# plot data
nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = treatment, colour = replicate), size = 5) +
  labs(title = "NMDS Species Percent",
       x = "NMDS1", 
       y = "NMDS2", 
       shape = "Treatment", 
       colour = "Replicate")

nmds_plot

# save plot as svg
# svg("NMDS_Species_Percent.svg")
# nmds_plot
# dev.off()
