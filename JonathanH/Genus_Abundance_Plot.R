# This script plots a stacked bar graph for samples containing percentage data.
# The input is genus clade percent data from Pavian (kraken reports)
# that have Eukaryota filtered out.

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
clean_data <- select(data, -taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")

# scale data and convert to data frame
scaled_data <- apply(clean_data[, -1],
                         MARGIN = 1,
                         FUN = function(x) x / sum(x)) %>%
  t() %>%
  as.data.frame()

# select for core taxa
core_data <- select(scaled_data,
                    "Gilliamella",
                    "Snodgrassella",
                    "Bifidobacterium",
                    "Lactobacillus",
                    "Frischella")

# TODO: setup treatment info
treatments <- rep(c("Control", "CLO", "THI"), 5)

# TODO: setup replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")

core_data$treatment <- treatments
core_data$replicate <- replicates

# TODO: for ordering on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
core_data$treatment <- factor(core_data$treatment,
                              levels = order)

# convert data frame into "long" format for stacked bar plot
long_data <- pivot_longer(core_data,
                          cols = 1:(ncol(core_data)-2),
                          names_to = "clade",
                          values_to = "percentage")

# plot data
# TODO: change title label
abundance_plot <- ggplot(long_data, aes(x = treatment,
                                        y = percentage,
                                        fill = clade)) +
  geom_bar(stat = "identity", colour = "black") +
  facet_grid(~replicate) +
  labs(title = "CTX Abundance Using Percent(%) Data",
       x = "Treatment",
       y = "Relative Abundance(%)",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1))

abundance_plot

# un-comment last 3 lines to save plot as svg
# TODO: change file name
# svg("CTX_Abundance_Plot.svg")
# abundance_plot
# dev.off()
