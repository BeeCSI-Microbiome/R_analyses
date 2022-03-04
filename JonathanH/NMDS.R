# Script for quick NMDS using the raw_clade table from R workflow
# This script is identical to what the R workflow would produce
# The purpose for this is to quickly test iterations without needing to 
# continuously run the workflow itself

# Script Author(s): Jonathan Ho
# Updated: Feb 4, 2022

library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(glue)


# User Defined Variables --------------------------------------------------
dataset_name <- 'clo_2020'
datapath <- 'results/clo_2020/plot_data/clo_raw_clade.csv'
treat_names <- c("Control", "Acute", "Sublethal")
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 4", "Rep 5")


# Functions ---------------------------------------------------------------
# cleans data into tidy format
tidy_data <- function(d) {
  clean_data <- select(d, -taxRank, -lineage) %>%
    pivot_longer(!name, names_to = "sample", values_to = "value") %>%
    pivot_wider(names_from = "name", values_from = "value")
  
  return(clean_data)
}

# calculates proportions and convert to data frame
calc_prop <- function(d) {
  sample_col <- select(d, "sample")
  prop_data <- select(d, -"sample") %>%
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample)
  
  return(prop_data)
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

# calculates nmds using Bray Curtis dissimilarity using only taxa data
calc_nmds <- function(data) {
  set.seed(1)
  sample_col <- select(data, "sample")
  
  nmds_data <- select(data, -"sample") %>%
    as.matrix() %>%
    metaMDS(distance = "bray") %>%
    scores() %>%
    as.data.frame()%>%
    mutate(sample_col) %>%
    relocate(sample)
  
  return(nmds_data)
}

# runs anosim and saves results in a text file in results folder
calc_ano <- function(d, group_data, taxa_level, dataset_name) {
  heading <- paste(taxa_level, "ANOSIM results:\n")
  ano_data <- select(d, -"sample") %>%
    as.matrix()
  
  rep_ano <- anosim(ano_data,
                    group_data$replicate,
                    distance = "bray",
                    permutations = 9999)
  
  treat_ano <- anosim(ano_data,
                      group_data$treatment,
                      distance = "bray",
                      permutations = 9999)
  
  cat(heading, file = glue("results/{dataset_name}/test_anosim.txt"), append = T)
  utils::capture.output(rep_ano, file = glue("results/{dataset_name}/test_anosim.txt"), append = T)
  utils::capture.output(treat_ano, file = glue("results/{dataset_name}/test_anosim.txt"), append = T)
}

# preps nmds data and runs ANOSIM on data
prep_and_ano <- function(d, treat_names, rep_names, taxa_level, dataset_name) {
  prop_data <- tidy_data(d) %>%
    calc_prop()
  
  group_data <- prop_data %>%
    calc_nmds() %>%
    treat_reps(treat_names, rep_names)
  
  calc_ano(prop_data, group_data, taxa_level, dataset_name)
  
  return(group_data)
}

# plots the nmds for activity 1 data
plot_nmds_1 <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  
  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(shape = 21, size = 3) +
    geom_polygon(data = hull, alpha = 0.5) +
    aes(fill = !!hull_var) +
    labs(title = plot_title,
         x = "NMDS1", 
         y = "NMDS2",
         fill = tools::toTitleCase(h_var))
  
  plot
}

# plots the nmds for activity 2 data
plot_nmds_2 <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  
  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data = hull, alpha = 0.5) +
    geom_point(aes(shape = replicate), size = 3) +
    aes(fill = !!hull_var) +
    labs(title = plot_title,
         x = "NMDS1", 
         y = "NMDS2",
         shape = "Replicate",
         fill = tools::toTitleCase(h_var))
  
  plot
}

# plots and saves nmds' with separate hulls for replicate and treatment
act_1_nmds <- function(genus_data, speci_data, dataset_name) {
  genus_treat_nmds <- plot_nmds_1(genus_data,
                                  "treatment",
                                  "Genus NMDS - Brays Curtis")
  genus_reps_nmds <- plot_nmds_1(genus_data,
                                 "replicate",
                                 "Genus NMDS - Brays Curtis")
  speci_treat_nmds <- plot_nmds_1(speci_data,
                                  "treatment",
                                  "Species NMDS - Bray Curtis")
  speci_reps_nmds <- plot_nmds_1(speci_data,
                                 "replicate",
                                 "Species NMDS - Bray Curtis")
  
  ggsave(plot = genus_treat_nmds, 
         filename = glue("results/{dataset_name}/test_genus_treat_nmds.png"), bg = "white")
  ggsave(plot = genus_reps_nmds, 
         filename = glue("results/{dataset_name}/test_genus_reps_nmds.png"), bg = "white")
  ggsave(plot = speci_treat_nmds, 
         filename = glue("results/{dataset_name}/test_speci_treat_nmds.png"), bg = "white")
  ggsave(plot = speci_reps_nmds, 
         filename = glue("results/{dataset_name}/test_speci_reps_nmds.png"), bg = "white")
}

# plots and saves nmds' for treatment only
act_2_nmds <- function(genus_data, speci_data, dataset_name) {
  genus_treat_nmds <- plot_nmds_2(genus_data,
                                  "treatment",
                                  "Genus NMDS - Brays Curtis")
  speci_treat_nmds <- plot_nmds_2(speci_data,
                                  "treatment",
                                  "Species NMDS - Bray Curtis")
  ggsave(plot = genus_treat_nmds,
         filename = glue("results/{dataset_name}/test_genus_treat_nmds.png"), bg = "white")
  ggsave(plot = speci_treat_nmds,
         filename = glue("results/{dataset_name}/test_speci_treat_nmds.png"), bg = "white")
}

# make and save nmds plots for genus and species levels using read data
# uses raw reads to calculate proportions
# export("make_nmds_plots")
make_nmds_plots <- function(data, treat_names, rep_names, dataset_name) {
  genus_data <- filter(data, taxRank == "G") %>%
    prep_and_ano(treat_names, rep_names, "Genus", dataset_name)
  speci_data <- filter(data, taxRank == "S") %>%
    prep_and_ano(treat_names, rep_names, "Species", dataset_name)
  
  utils::write.csv(genus_data,
                   file = glue("results/{dataset_name}/plot_data/test_genus_nmds.csv"),
                   row.names = F)
  utils::write.csv(speci_data,
                   file = glue("results/{dataset_name}/plot_data/test_species_nmds.csv"),
                   row.names = F)
  
  # split act 1 and 2 nmds here, based on number of treatments
  if (length(treat_names) > 2) {
    # can make convex hulls for both replicate and treatment
    act_1_nmds(genus_data, speci_data, dataset_name)
  } else {
    # can only make hulls for treatment
    act_2_nmds(genus_data, speci_data, dataset_name)
  }
}


# Call --------------------------------------------------------------------
data <- read_csv(datapath)
make_nmds_plots(data, treat_names, rep_names, dataset_name)
