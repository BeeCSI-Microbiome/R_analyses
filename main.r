# Author(s): Lance Lansing, Jonathan Ho
# Nov 2021
# Main script for analysis of taxonomic classification results
#
# Input: table from Pavian with clades AND taxon counts, not collapsed

setwd("~/beecsi/R_analyses")

# ------------------------------- Package Setup --------------------------------
packages <- c('tidyverse', 'vegan', 'modules', "data.table", "metagenomeSeq",
              "ggplot2", "glue")
lapply(packages, library, character.only = TRUE)
# ______________________________________________________________________________

# ------------------------- Load Aux Files as Modules --------------------------
rsummary <- use("scripts/reads_summary.R") 
ip <- use('scripts/initial_processing.R')
scaling <- use('scripts/reads_scaling.R')
exploratory <- use('scripts/exploratory_functions.R')
# ______________________________________________________________________________


# ---------------------------------- Globals -----------------------------------
# dataset name for file writing purposes
dataset_name <- "cor_2020"
# create result folders using the dataset name
ifelse(!dir.exists(glue("results/{dataset_name}")),
       dir.create(glue("results/{dataset_name}"), mode='777'), FALSE)
ifelse(!dir.exists(glue("results/{dataset_name}/plot_data")),
       dir.create(glue("results/{dataset_name}/plot_data"), mode='777'), FALSE)

# Input file paths:
# Counts table (clade and taxon counts, uncollapsed)
counts_path <- "../data/cor_2020/cor_all_taxa.tsv"
 
# Treatment and replicate names
<<<<<<< Updated upstream
treat_names <- c("Control", "CLO", "THI")
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 4", "Rep 5")

=======
treat_names <- c("exposed", "unexposed")
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 4", "Rep 5")
>>>>>>> Stashed changes
# Percentile value used by CSS (default=0.5)
css_percentile = 0.5
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
<<<<<<< Updated upstream
                                    rep_names)

# Do not use, only for CTX experiment
# exploratory$make_separate_ctx_bars(tables[["raw_clade"]],
#                                    treat_names,
#                                    rep_names)
=======
                                    rep_names,
                                    dataset_name)
>>>>>>> Stashed changes

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

