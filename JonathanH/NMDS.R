# This script produces an NMDS plot using percentage data.
# The input is percent data downloaded from the Comparison tab in Pavian
# without collapsing taxa.

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, un-collapse the taxa, and download.

# This script assumes that samples are grouped by replicates first, then treatments.
# Example below:
# Rep 1 TreatmentA, Rep 1 TreatmentB, Rep 1 TreatmentC, 
# Rep 2 TreatmentA, Rep 2 TreatmentB, Rep 2 TreatmentC, ...

# TODO:'s show recommended fields that should be changed for each analysis

# Script Author(s): Jonathan Ho

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)

# TODO: change file path
datapath <- '../2020_clo_kraken2/clo_kraken_all_rawread_uncollapsed.tsv'
# TODO: setup treatment info
treat_names <- c("Control", "Acute", "Sublethal")
# TODO: setup replicate info
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 5", "Rep 6")


# read data
data <- read.delim(file = datapath,
                       header = TRUE,
                       sep = '\t')

# clean data and filter out non-bacteria
clean_data <- filter(data, taxRank == "G") %>%
  select(-taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")%>%
  select(-"Apis",
         -"Aspergillus",
         -"Ascosphaera", 
         -"Nosema", 
         -"Crithidia", 
         -"Lotmaria")

# scale data and convert to data frame
scaled_data <- apply(clean_data[,2:ncol(clean_data)],
                     MARGIN = 1,
                     FUN = function(x) x / sum(x) * 100) %>%
  t() %>%
  as.data.frame()

# convert taxa part of data frame to matrix for nmds
taxa_mat <- as.matrix(scaled_data)

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

# nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) + 
#   geom_point(aes(colour = replicate, shape = treatment), size = 5) +
#   labs(title = "NMDS - Bray Curtis",
#        x = "NMDS1", 
#        y = "NMDS2",
#        colour = "Replicate",
#        shape = "Treatment")

# calculate convex hull
hull <- nmds_data %>%
  group_by(treatment) %>%
  slice(chull(NMDS1, NMDS2))

# make nmds plot
nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(shape = 21, size = 3) +
  aes(fill = treatment) +
  geom_polygon(data = hull, alpha = 0.5) +
  labs(title = "NMDS - Bray Curtis",
       x = "NMDS1", 
       y = "NMDS2")

nmds_plot




# un-comment last 3 lines to save plot as svg
# svg("NMDS_Percent.svg")
# nmds_plot
# dev.off()


################################################

# stats: ANOSIM

ano = anosim(taxa_mat,
             nmds_data$replicate,
             distance = "bray",
             permutations = 9999)

ano
