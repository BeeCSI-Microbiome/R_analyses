# This script plots a bar graph for samples containing percentage data.
# The input is genus clade percent data from Pavian (kraken reports).
# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, filter out Eukaryota, ensure clade is selected and select percent.

library(dplyr)
library(tidyr)
library(ggplot2)

# read data
ctx_data <- read.delim(file = 'ctx_kraken_genus_percent_data.tsv',
                       header = TRUE,
                       sep = '\t')

# clean data
bifido_data <- select(ctx_data, -taxRank, -taxID, -Max, -lineage) %>%
  filter(name == "Bifidobacterium") %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")

# add treatment info
treatments <- c("Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI")
bifido_data$treatment <- treatments

# for ordering on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
bifido_data$treatment <- factor(bifido_data$treatment,
                                   levels = order)

# add replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")
bifido_data$replicate <- replicates

bifido_plot <- ggplot(bifido_data, 
                      aes(x = replicate, 
                          y = Bifidobacterium, 
                          fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(limits = order) +
  labs(title = "Bifidobacterium Percent Abundance",
       x = "Replicate",
       y = "Percent",
       fill = "Treatment")

bifido_plot

# save plot as svg
# svg("Bifido_Percent_Abundance.svg")
# bifido_plot
# dev.off()
