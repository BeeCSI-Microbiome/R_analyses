# This script plots a bar graph for a single genus in samples containing percentage data.
# The input is genus clade percent data from Pavian (kraken reports).

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, filter out Eukaryota, ensure clade is selected and select percent.

# TODO:'s show recommended fields that should be changed for each analysis

library(dplyr)
library(tidyr)
library(ggplot2)

# TODO: change file path
datapath <- '2020_ctx_kraken2/ctx_kraken_genus_percent.tsv'

# read data
data <- read.delim(file = datapath,
                       header = TRUE,
                       sep = '\t')

# clean data
# TODO: change name to be genus of choice
data <- select(data, -taxRank, -taxID, -Max, -lineage) %>%
  filter(name == "Bifidobacterium") %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")

# TODO: add treatment info
treatments <- rep(c("Control", "CLO", "THI"), 5)
data$treatment <- treatments

# TODO: for ordering on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
data$treatment <- factor(data$treatment,
                                   levels = order)

# TODO: add replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")
data$replicate <- replicates

# plot data
# TODO: change x to be genus of choice and title label
genus_plot <- ggplot(data,
                     aes(x = replicate,
                         y = Bifidobacterium,
                         fill = treatment))+
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(limits = order) +
  labs(title = "Bifidobacterium Percent Abundance",
       x = "Replicate",
       y = "Percent",
       fill = "Treatment")

genus_plot

# un-comment last 3 lines to save plot as svg
# TODO: change file name
# svg("Bifido_Percent_Abundance.svg")
# genus_plot
# dev.off()
