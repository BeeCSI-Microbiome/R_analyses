# Script for running ANCOM-BC for for pairwise differential 
# abundance analysis and producing volcano plots.
# Assumptions: using treatment+replicate for formula
# Takes two functions from exploratory_functions.R for data wrangling:
# tidy_data and treat_reps
# Takes as input: raw clade counts from entire data set

# Author(s): Jonathan Ho
# Updated: Apr 11, 2022

library(tidyverse)
library(ANCOMBC)
library(microbiome)

# User Defined Variables --------------------------------------------------
datapath <- 'results/thi_2020/plot_data/thi_raw_clade.csv'
treat_names <- c("Control","Acute", "Sublethal")
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")
plot_title <- "Caged Control vs Acute THI - Genus DA"
plot_file_name <- 'thi_control_acute_genus_treatrep.png'
summary_file_name <- "thi_2020_genus_mod2_trunc.csv"


# Functions from exploratory_functions.R ----------------------------------

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
# do both genus and species level
data <- read_csv(datapath) %>%
  filter(taxRank == "G")

abun_data <- select(data, -taxRank, -lineage) %>%
  column_to_rownames('name') %>%
  as.matrix()

metadata <- tidy_data(data) %>%
  select(sample) %>%
  treat_reps(treat_names, rep_names) %>%
  column_to_rownames('sample') %>%
  add_labs(lab_names, treat_names)

OTU <- otu_table(abun_data, taxa_are_rows = T)
samples <- sample_data(metadata)

phylo_obj <- phyloseq(OTU, samples)


