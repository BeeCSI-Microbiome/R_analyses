# This script produces an NMDS plot using species percentage data.
# The input is species percent data from Pavian (kraken reports).

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, filter out Eukaryota, ensure clade is selected and select percent.

# TODO:'s show recommended fields that should be changed for each analysis

# Leave in more info? bees and eukaryotes to capture more info?

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# TODO: change file path
datapath <- '2020_ctx_kraken2/ctx_kraken_genus_percent.tsv'

# read data
data <- read.delim(file = datapath,
                       header = TRUE,
                       sep = '\t')

# clean data
clean_data <- select(data, -taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")

# convert taxa part of data frame to matrix for nmds
taxa_data <- clean_data[,2:ncol(clean_data)]
taxa_mat <- as.matrix(taxa_data) 

# perform nmds
set.seed(1)
nmds = metaMDS(taxa_mat, distance = "bray")

# TODO: setup treatment info
treatments <- rep(c("Control", "CLO", "THI"), 5)

# TODO: setup replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")

# add back information to matrix
nmds_data <- as.data.frame(scores(nmds))
nmds_data$sample <- clean_data$sample
nmds_data$replicate <- replicates
nmds_data$treatment <- treatments

# TODO: for ordering treatments on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
nmds_data$treatment <- factor(nmds_data$treatment,
                              levels = order)

# plot data
# TODO: change title label
nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = treatment, colour = replicate), size = 5) +
  labs(title = "CTX NMDS Species Percent",
       x = "NMDS1", 
       y = "NMDS2", 
       shape = "Treatment", 
       colour = "Replicate")

nmds_plot

# un-comment last 3 lines to save plot as svg
# TODO: change file name
# svg("CTX_NMDS_Species_Percent.svg")
# nmds_plot
# dev.off()
