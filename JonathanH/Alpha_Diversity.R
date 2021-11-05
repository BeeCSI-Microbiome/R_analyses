# This script calculates and plots alpha diversity metrics for 
# kraken output data downloaded through Pavian.
# The input is raw read data downloaded from the Comparison tab
# without collapsing taxa.

# To download the data in a script-ready form, go to the Comparison tab in
# Pavian, un-collapse the taxa, and download.

# This script assumes that samples are grouped by replicates first, then treatments.
# Example below:
# Rep 1 TreatmentA, Rep 1 TreatmentB, Rep 1 TreatmentC, Rep 2 TreatmentA,
# Rep 2 TreatmentB, Rep 2 TreatmentC, ...

# TODO:'s show recommended values that should be changed for each analysis

# Script Author(s): Jonathan Ho, Jocelyn O'brien

library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)

# TODO: change file path
datapath <- '2020_clo_kraken2/clo_kraken_all_rawread_uncollapsed.tsv'

# TODO: setup treatment info
treat_names <- c("Control", "Acute", "Sublethal")

# TODO: setup replicate info
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 5", "Rep 6")

# read data
data <- read.delim(file = datapath,
                   header = TRUE,
                   sep = '\t')

# filter for species and clean data
clean_data <- filter(data, taxRank == "S") %>%
  select(-taxRank, -taxID, -Max, -lineage) %>%
  pivot_longer(!name, names_to = "sample", values_to = "percent") %>%
  pivot_wider(names_from = "name", values_from = "percent") %>%
  select(-"Apis mellifera")

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

alpha_df <- calc_diversity_df(clean_data[,2:ncol(clean_data)])

# add in treatment and replicate cols
num_treats <- length(treat_names)
num_reps <- length(rep_names)

treatments <- rep(treat_names, num_reps)
replicates <- c()
for (r in 1:num_reps) {
  replicates = c(replicates, rep(rep_names[r], num_treats))
}

alpha_df$treatment <- treatments
alpha_df$replicate <- replicates

# adjust factor levels for ordering
alpha_df$treatment <- factor(alpha_df$treatment,
                                levels = treat_names)

# make alpha diversity plot
alpha_plot <- ggplot(alpha_df, aes(x = treatment, y = Inv_Simpson)) +
  geom_boxplot() +
  labs(title = "Alpha Diversity Plot",
       x = "Treatment")

alpha_plot

# un-comment last 3 lines to save plot as svg
# svg("Alpha_Diversity_Plot.svg")
# abundance_plot
# dev.off()

# alpha_df$sample <- gsub('(\\w+\\d\\d(u|e)).*$', "\\1", acp$sample)

######################################################

# Stats: Mann Whitney U Test

# test between control and acute treatments
control_vs_acute <- filter(alpha_df, treatment %in% c("Control", "Acute"))
wilcox.test(Inv_Simpson~treatment, data=control_vs_acute)

# test between control and sublethal treatments
control_vs_sublethal <- filter(alpha_df, treatment %in% c("Control", "Sublethal"))
wilcox.test(Inv_Simpson~treatment, data=control_vs_sublethal)

# test between acute and sublethal treatments
acute_vs_sublethal <- filter(alpha_df, treatment %in% c("Acute", "Sublethal"))
wilcox.test(Inv_Simpson~treatment, data=acute_vs_sublethal)
