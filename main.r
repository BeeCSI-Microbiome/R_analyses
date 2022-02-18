# Author(s): Lance Lansing, Jonathan Ho
# Nov 2021
# Main script for analysis of taxonomic classification results
#
# Input: table from Pavian with clades AND taxon counts, not collapsed

setwd("~/beecsi/R_analyses")

# ------------------------------- Package Setup --------------------------------
# "here", "PMCMRplus", "broom", "statmod" added for differential_abundance.R <DA>
packages <- c('tidyverse', 'vegan', 'modules', "data.table", "metagenomeSeq",
              "ggplot2", "glue", "here", "PMCMRplus", "broom", "statmod")
lapply(packages, library, character.only = TRUE)

# Get utility functions <DA> ---------------------------------------------------
source(here::here('scripts','meg_utility_functions.R'))

# ______________________________________________________________________________

# ------------------------- Load Aux Files as Modules --------------------------
rsummary <- use("scripts/reads_summary.R") 
ip <- use('scripts/initial_processing.R')
scaling <- use('scripts/reads_scaling.R')
exploratory <- use('scripts/exploratory_functions.R')
widen_results <- use('scripts/widenResults.R')
differential_abundance_kraken <- use('scripts/differential_abundance.R')
# ______________________________________________________________________________


# ---------------------------------- Globals -----------------------------------
# dataset name for file writing purposes
dataset_name <- "cor_2020"
# create result folders using the dataset name
ifelse(!dir.exists(glue("results/{dataset_name}")),
       dir.create(glue("results/{dataset_name}"), mode='777'), FALSE)
ifelse(!dir.exists(glue("results/{dataset_name}/plot_data")),
       dir.create(glue("results/{dataset_name}/plot_data"), mode='777'), FALSE)
# <DA> add sub directory for DA analysis
ifelse(!dir.exists(glue("results/{dataset_name}/aggregated_dir")),
       dir.create(glue("results/{dataset_name}/aggregated_dir"), mode='777'), FALSE)
# <DA> aggreagated data directory for use by differential_abundance.R
# aggregated_dir <- here("aggregated_data_for_analysis") # produced by widen results should be in func ARG. 

ifelse(
  !dir.exists(aggregated_dir), 
  dir.create((aggregated_dir), mode='777'), 
  FALSE
)

# Input file paths:
# Counts table (clade and taxon counts, uncollapsed)
counts_path <- "../data/cor_2020/cor_all_taxa.tsv"
metadata_filepath = "../data/cor_2020/cor_2020_metadata.csv" # can this be better? <DA>

# Kraken Paths <DA>
krakenReportPaths <- Sys.glob("../data/hbb_2020/kraken_reports/*_report.txt") # set to global
krakenReportNames <- list.files(path = "../data/_2020/kraken_reports/", pattern = '*_report.txt')

# <DA> unsure if this is correct -kurt. 
kraken_analytical <- Sys.glob(here::here("aggregated_data_for_analysis", dataset_name, "krakenAnalytical_*.csv"))

# Metadata
metadata <- read.csv(metadata_filepath, header=T)
 
# Treatment and replicate names
treat_names <- c("exposed", "unexposed")
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 4", "Rep 5")

# Percentile value used by CSS (default=0.5)
css_percentile = 0.5

# Output paths <DA>; can these just be placed in internal functions? 
# Set the output directory for statistics:
stats_output_dir = 'results/differential_abundance_stats'

# Create output directories if they don't exist
ifelse(!dir.exists(graph_output_dir), dir.create((graph_output_dir), mode='777'), FALSE)
ifelse(!dir.exists(stats_output_dir), dir.create((stats_output_dir), mode='777'), FALSE)
# ______________________________________________________________________________


# -------------------------------- Read Input ----------------------------------
# Read in count table
ct <- read_tsv(counts_path)
# ______________________________________________________________________________


# ------------------------------- Reads Summary --------------------------------
reads_summary <- rsummary$create_summary_table(ct, dataset_name)
# ______________________________________________________________________________


# -------------------- Formatting, Filtering, and Scaling ----------------------
# Format and perform filtering on the table
ct <- ip$format_count_table(ct) %>%
  ip$filter_table() %>%
  ip$group_taxa_of_interest()

tables <- scaling$scaling_procedure(ct, css_percentile, dataset_name)

mrexp <- tables$css_MRexp
raw_taxon <- tables$raw_taxon
raw_clade <- tables$raw_clade
scaled_taxon <- tables$scaled_taxon
scaled_clade <- tables$scaled_clade
# ______________________________________________________________________________



# Relative Abundance ------------------------------------------------------
exploratory$make_interest_abundance(tables[["raw_clade"]],
                                    treat_names,
                                    rep_names,
                                    dataset_name,
                                    additional_taxa)


# Alpha Diversity ---------------------------------------------------------
exploratory$make_all_alpha_plots(tables[["raw_clade"]],
                                 treat_names,
                                 rep_names,
                                 dataset_name)


# Beta Diversity ----------------------------------------------------------
exploratory$make_nmds_plots(tables[["raw_clade"]],
                            treat_names,
                            rep_names,
                            dataset_name)


# Differential Abundance ----------------------------------------------------------
widen_results$widen_results_function(krakenReportPaths)
differential_abundance_kraken$kraken_differential_abundance(dataset_name)
# ______________________________________________________________________________







