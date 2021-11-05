# This script plots a stacked bar graph for samples containing percentage data.
# The input is percent data downloaded from the Comparison tab in Pavian
# without collapsing taxa.

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, un-collapse the taxa, and download.

# This script assumes that samples are grouped by replicates first, then treatments.
# Example below:
# Rep 1 TreatmentA, Rep 1 TreatmentB, Rep 1 TreatmentC, 
# Rep 2 TreatmentA, Rep 2 TreatmentB, Rep 2 TreatmentC, ...

# TODO:'s show recommended values that should be changed for each analysis

# Script Author(s): Jonathan Ho, Gale Chen

library(dplyr)
library(tidyr)
library(ggplot2)

# TODO: change file path
datapath <- '2020_ctx_kraken2/ctx_kraken_all_percent_uncollapsed.tsv'

# TODO: setup treatment info
treat_names <- c("Control", "CLO", "THI")

# TODO: setup replicate info
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")

# TODO: give a title for the plot
plot_title <- "CTX Abundance Using Percent(%) Data"


# read data
data <- read.delim(file = datapath,
                       header = TRUE,
                       sep = '\t')

# filter to genus level and clean data
clean_data <- filter(data, taxRank == "G") %>%
  select(-taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")%>%
  select(-Apis)

# scale data and convert to data frame
scaled_data <- apply(clean_data[,2:ncol(clean_data)],
                     MARGIN = 1,
                     FUN = function(x) x / sum(x)) %>%
  t() %>%
  as.data.frame()

# select for core taxa and add Others col
core_data <- select(scaled_data,
                    "Gilliamella",
                    "Snodgrassella",
                    "Bifidobacterium",
                    "Lactobacillus",
                    "Frischella") %>%
  mutate(Others = 1 - (Gilliamella +
                         Snodgrassella + 
                         Bifidobacterium +
                         Lactobacillus +
                         Frischella))

# add in treatment and replicate cols
num_treats <- length(treat_names)
num_reps <- length(rep_names)

treatments <- rep(treat_names, num_reps)
replicates <- c()
for (r in 1:num_reps) {
  replicates = c(replicates, rep(rep_names[r], num_treats))
}

core_data$treatment <- treatments
core_data$replicate <- replicates

# adjust factor levels for ordering
core_data$treatment <- factor(core_data$treatment,
                              levels = treat_names)


# convert data frame into "long" format for stacked bar plot
long_data <- pivot_longer(core_data,
                          cols = 1:(ncol(core_data)-2),
                          names_to = "clade",
                          values_to = "percentage")

# plot data
abundance_plot <- ggplot(long_data, 
                         aes(x = treatment,
                             y = percentage,
                             fill = factor(clade,
                                           levels = c("Bifidobacterium",
                                                      "Frischella",
                                                      "Gilliamella",
                                                      "Lactobacillus",
                                                      "Snodgrassella",
                                                      "Others")))) +
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
# svg("CTX_Abundance_Plot.svg")
# abundance_plot
# dev.off()
