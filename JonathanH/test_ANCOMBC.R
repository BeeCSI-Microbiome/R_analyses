# This script is for testing ANCOM_BC as another method for differential
# abundance. This is to be done before adding to the main R workflow.
# Takes some functions from exploratory_functions.R

# Author(s): Jonathan Ho
# Updated: Mar 19, 2022

library(tidyverse)
library(ANCOMBC)
library(microbiome)

# User Defined Variables --------------------------------------------------
dataset_name <- 'ctx_2020'
datapath <- 'results/ctx_2020/plot_data/ctx_raw_clade.csv'
treat_names <- c("Control", "CLO")
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")


# Functions ---------------------------------------------------------------
# cleans data into tidy format
tidy_data <- function(d) {
  clean_data <- select(d, -taxRank, -lineage) %>%
    pivot_longer(!name, names_to = "sample", values_to = "value") %>%
    pivot_wider(names_from = "name", values_from = "value")
  
  return(clean_data)
}

# add in treatment and replicate cols
# assumes same number of replicates for each treatment and vice versa. 
# assumes that samples are grouped by replicates first, then treatments.
# Example below:
# Rep 1 TreatmentA, Rep 1 TreatmentB, Rep 1 TreatmentC, 
# Rep 2 TreatmentA, Rep 2 TreatmentB, Rep 2 TreatmentC, ...
# Can double check that these are correct by comparing with samples col
treat_reps <- function(d, treat_names, rep_names) {
  num_treats <- length(treat_names)
  num_reps <- length(rep_names)
  
  treatments <- rep(treat_names, num_reps)
  replicates <- c()
  for (r in 1:num_reps) {
    replicates = c(replicates, rep(rep_names[r], num_treats))
  }
  
  d$treatment <- treatments
  d$replicate <- replicates
  
  # adjust factor levels for future plotting
  d$treatment <- factor(d$treatment,
                        levels = treat_names)
  return(d)
}


# Setup -------------------------------------------------------------------

data <- read_csv(datapath) %>%
  select(-contains('_dT')) %>%
  filter(taxRank == "G")

abun_data <- select(data, -taxRank, -lineage) %>%
  column_to_rownames('name') %>%
  as.matrix()

metadata <- tidy_data(data) %>%
  select(sample) %>%
  treat_reps(treat_names, rep_names) %>%
  column_to_rownames('sample')

OTU <- otu_table(abun_data, taxa_are_rows = T)
samples <- sample_data(metadata)

phylo_obj <- phyloseq(OTU, samples)


# Call --------------------------------------------------------------------
# 2 different models differing in their formulas
mod1 <- ancombc(phylo_obj,
                formula = 'treatment',
                p_adj_method = "BH",
                group = 'treatment',
                global = T)

# includes replicate as confounding variable
mod2 <- ancombc(phylo_obj,
                formula = 'treatment+replicate',
                p_adj_method = "BH",
                group = 'treatment',
                global = T)


write.csv(mod1$res,
          file = "ctx_2020_clo_genus_mod1.csv",
          row.names = T)

write.csv(mod2$res,
          file = "ctx_2020_clo_genus_mod2.csv",
          row.names = T)


mod1$res
mod1$res_global

mod2$res
mod2$res_global
