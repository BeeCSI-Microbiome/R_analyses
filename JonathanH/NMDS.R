# This script produces an NMDS plot using percentage data.
# The input is percent data downloaded from the Comparison tab in Pavian
# without collapsing taxa.

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, un-collapse the taxa, and download.

# TODO:'s show recommended fields that should be changed for each analysis

# Script Author(s): Jonathan Ho

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# TODO: change file path
datapath <- '2020_ctx_kraken2/ctx_kraken_all_percent_uncollapsed.tsv'

# TODO: setup treatment info
treat_names <- c("Control", "CLO", "THI")

# TODO: setup replicate info
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")

# TODO: give a title for the plot
plot_title <- "CTX Bifidobacterium Percent Abundance"

# read data
data <- read.delim(file = datapath,
                       header = TRUE,
                       sep = '\t')

# clean data
clean_data <- filter(data, taxRank == "G") %>%
  select(-taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")%>%
  select(-Apis)

# convert taxa part of data frame to matrix for nmds
taxa_data <- clean_data[,2:ncol(clean_data)]
taxa_mat <- as.matrix(taxa_data) 

# perform nmds
set.seed(1)
nmds = metaMDS(taxa_mat, distance = "bray")

# add back information to matrix
nmds_data <- as.data.frame(scores(nmds))
nmds_data$sample <- clean_data$sample

num_treats <- length(treat_names)
num_reps <- length(rep_names)

treatments <- rep(treat_names, num_reps)
replicates <- c()
for (r in 1:num_reps) {
  replicates = c(replicates, rep(rep_names[r], num_treats))
}

nmds_data$treatment <- treatments
nmds_data$replicate <- replicates

# adjust factor levels for ordering
nmds_data$treatment <- factor(nmds_data$treatment,
                              levels = treat_names)

# plot data
nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(shape = treatment, colour = replicate), size = 5) +
  labs(title = plot_title,
       x = "NMDS1", 
       y = "NMDS2", 
       shape = "Treatment", 
       colour = "Replicate")

nmds_plot

# un-comment last 3 lines to save plot as svg
# svg("CTX_NMDS_Species_Percent.svg")
# nmds_plot
# dev.off()
