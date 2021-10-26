# This script plots a stacked bar graph for samples containing percentage data.
# The input is genus clade percent data from Pavian (kraken reports)
# that have Eukaryota filtered out.

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, filter out Eukaryota, ensure clade is selected and select percent.

# TODO:'s show recommended values that should be changed for each analysis

# This script assumes that samples are grouped by replicates first, then treatments.
# Example below:
# Rep 1 TreatmentA, Rep 1 TreatmentB, Rep 1 TreatmentC, Rep 2 TreatmentA,
# Rep 2 TreatmentB, Rep 2 TreatmentC, ...

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

# TODO: change file path
datapath <- '2020_clo_kraken2/clo_kraken_species_rawread.tsv'

# TODO: setup treatment info
treat_names <- c("Control", "Acute", "Sublethal")

# TODO: setup replicate info
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 5", "Rep 6")

# read data
data <- read.delim(file = datapath,
                   header = TRUE,
                   sep = '\t')

# clean data
clean_data <- select(data, -taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent")

# calc alpha metrics
calc_diversity_df <- function(x){
  observed_richness <- specnumber(x)
  invsimpson <- diversity(x, index="invsimpson")
  simpson <- diversity(x, index="simpson")
  shannon <- diversity(x, index="shannon")
  evenness <- shannon/log(observed_richness)
  div_df <- data.frame(
    ID = clean_data$sample,
    Observed_Richness = observed_richness,
    Inv_Simpson = invsimpson,
    Simpson = simpson,
    Shannon = shannon,
    Evenness = evenness
  )
  div_df
}

alpha_div <- calc_diversity_df(clean_data[,2:ncol(clean_data)])

# add in treatment and replicate cols
num_treats <- length(treat_names)
num_reps <- length(rep_names)

treatments <- rep(treat_names, num_reps)
replicates <- c()
for (r in 1:num_reps) {
  replicates = c(replicates, rep(rep_names[r], num_treats))
}

alpha_div$treatment <- treatments
alpha_div$replicate <- replicates

# adjust factor levels for ordering
alpha_div$treatment <- factor(alpha_div$treatment,
                                levels = treat_names)

# acp$sample <- gsub('(\\w+\\d\\d(u|e)).*$', "\\1", acp$sample)

alpha_plot <- ggplot(alpha_div, aes(x = treatment, y = Shannon)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity Plot",
       x = "Treatment",
       y = "Index")

alpha_plot

