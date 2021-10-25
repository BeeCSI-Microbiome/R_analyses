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
