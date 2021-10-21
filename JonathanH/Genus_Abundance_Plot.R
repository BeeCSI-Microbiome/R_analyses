# This script plots a stacked bar graph for samples containing percentage data.
# The input is genus clade percent data from Pavian (kraken reports)
# that have Eukaryota filtered out.

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, filter out Eukaryota, ensure clade is selected and select percent.

# TODO:'s show recommended values that should be changed for each analysis

# This script assumes that samples are grouped by replicates first, then treatments.
# example below:
# Rep 1 TreatmentA
# Rep 1 TreatmentB
# Rep 1 TreatmentC
# Rep 2 TreatmentA
# Rep 2 TreatmentB
# Rep 2 TreatmentC

library(dplyr)
library(tidyr)
library(ggplot2)

# TODO: change file path
datapath <- '2020_ctx_kraken2/ctx_kraken_genus_percent.tsv'

# TODO: setup treatment info
treatment_vec <- c("Control", "CLO", "THI")
# TODO: number of replicates
num_reps <- 5

# TODO: setup replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")

# TODO: Give a title for the plot
plot_title <- "CTX Abundance Using Percent(%) Data"

# read data
data <- read.delim(file = datapath,
                       header = TRUE,
                       sep = '\t')

# clean data
clean_data <- select(data, -taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")

# select for core taxa
core_data <- select(clean_data,
                    "Gilliamella",
                    "Snodgrassella",
                    "Bifidobacterium",
                    "Lactobacillus",
                    "Frischella")

# scale data and convert to data frame
scaled_data <- apply(core_data,
                     MARGIN = 1,
                     FUN = function(x) x / sum(x)) %>%
  t() %>%
  as.data.frame()

# add in treatment and replicate cols
treatments <- rep(treatment_vec, num_reps)
scaled_data$treatment <- treatments
scaled_data$replicate <- replicates

# adjust factor levels for ordering
scaled_data$treatment <- factor(scaled_data$treatment,
                              levels = treatment_vec)

# convert data frame into "long" format for stacked bar plot
long_data <- pivot_longer(scaled_data,
                          cols = 1:(ncol(scaled_data)-2),
                          names_to = "clade",
                          values_to = "percentage")

# plot data
abundance_plot <- ggplot(long_data, aes(x = treatment,
                                        y = percentage,
                                        fill = clade)) +
  geom_bar(stat = "identity", colour = "black") +
  facet_grid(~replicate) +
  labs(title = plot_title,
       x = "Treatment",
       y = "Relative Abundance(%)",
       fill = "Genus") +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1))

abundance_plot

# un-comment last 3 lines to save plot as svg
# svg("Genus_Abundance_Plot.svg")
# abundance_plot
# dev.off()
