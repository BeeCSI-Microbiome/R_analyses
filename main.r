# Author(s): Lance Lansing, Jonathan Ho, Kurt Clarke

# Main script for analysis of taxonomic classification results

# Input: table from Pavian with clades AND taxon counts, not collapsed
# User defined inputs (eg input file paths) must be set under section "Globals"

setwd("~/beecsi/R_analyses")

# Package setup -----------------------------------------------------------
packages <- c("tidyverse",
              "vegan",
              "modules",
              "data.table",
              "metagenomeSeq",
              "ggplot2",
              "glue")
lapply(packages, library, character.only = TRUE)
# _________________________________________________________________________


# Load aux scripts as modules ---------------------------------------------
rsummary <- use("scripts/reads_summary.R")
ip <- use("scripts/initial_processing.R")
exploratory <- use("scripts/exploratory_functions.R")
da_ancombc <- use("scripts/da-ancombc.R")
indicsp <- use("scripts/indicator_taxa_analysis.R")
# _________________________________________________________________________


# Globals -----------------------------------------------------------------
# Name of the dataset for file writing purposes
dataset_name <- "oxy_2021"

# Create strings for output directories
main_outdir <- glue("results/{dataset_name}")
nmds_dir <- glue("{main_outdir}/nmds_anosim")
alpha_div_dir <- glue("{main_outdir}/alpha_diversity")
rel_abund_dir <-  glue("{main_outdir}/relative_abundance")

# <DA> Subdirectory for DA results
da_dir <- glue("{main_outdir}/differential_abundance")
da_ancombc_dir <- glue("{da_dir}/ancombc")
ind_sp_dir <- glue("{da_dir}/indicator_species_analysis")
ind_sp_dir_clade <- glue("{ind_sp_dir}/raw_clade")
ind_sp_dir_taxon <- glue("{ind_sp_dir}/raw_taxon")

create_dir_if_nonexistant <- function(path) {
  ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
}
## Create output directories ####
lapply(c(main_outdir, nmds_dir, alpha_div_dir, rel_abund_dir, da_dir,
         da_ancombc_dir,ind_sp_dir, ind_sp_dir_clade, ind_sp_dir_taxon),
       create_dir_if_nonexistant)

## User-defined values ####
# Dataset-specific taxa of interest (provide a list of taxa strings)
# e.g. c("Lactobacillus", "Gilliamella apis")
additional_taxa <- c("Pantoea agglomerans")


# The following key should match the substring in the sample name that specifies 
# a treatment with the name of that treatment, including control
treatment_key <- list(d0="control",d1="oxytetracycline")
# a treatment with the name of that treatment
# THE CONTROL TREATMENT MUST BE THE FIRST ELEMENT (e.g. control, unexposed)
# _________________________________________________________________________


# Read input files --------------------------------------------------------
# Count table
ct <- read_csv(counts_path)
# _________________________________________________________________________


# Reads summary -----------------------------------------------------------
# Creates a table containing total, unclassified, and classified read counts,
# and % of reads classified as A. mellifera, Bacteria, and other specific taxa
reads_summary <- rsummary$create_summary_table(ct, dataset_name)
# _________________________________________________________________________


# # Differential Abundance --------------------------------------------------

# Formatting, filtering, and calculating clade counts ---------------------
# Format and perform filtering on the table
ct <- ip$format_count_table(ct) %>%
  ip$filter_table() %>%
  ip$group_taxa_of_interest()

tables <- ip$calculate_clade_counts(ct)

write.csv(tables[["raw_clade"]],
          glue("{main_outdir}/raw_clade_counts.csv"),
          row.names = FALSE)
write.csv(tables[["raw_taxon"]],
          glue("{main_outdir}/raw_taxon_counts.csv"),
          row.names = FALSE)

# source("scripts/snippets/taxa_cutoff_explore.R")
# __________________________________________________________________________

# ANCOMBC Differential Abundance ------------------------------------------
da_ancombc$run_ancombc(tables[["raw_clade"]],
                       treatment_key,
                       dataset_name,
                       "S",
                       da_ancombc_dir)
da_ancombc$run_ancombc(tables[["raw_clade"]],
                       treatment_key,
                       dataset_name,
                       "G",
                       da_ancombc_dir)

# Indicator Taxa Analysis--------------------------------------------------
indicsp$run_indicator_analysis(dataset_name,
                               tables[["raw_taxon"]],
                               ind_sp_dir_taxon,
                               treatment_key)
indicsp$run_indicator_analysis(dataset_name,
                               tables[["raw_clade"]],
                               ind_sp_dir_clade,
                               treatment_key)
indicsp$run_indicator_analysis(dataset_name,
                               tables[["raw_taxon"]],
                               ind_sp_dir_taxon,
                               treatment_key,
                               "S")
indicsp$run_indicator_analysis(dataset_name,
                               tables[["raw_clade"]],
                               ind_sp_dir_clade,
                               treatment_key,
                               "S")
indicsp$run_indicator_analysis(dataset_name,
                               tables[["raw_taxon"]],
                               ind_sp_dir_taxon,
                               treatment_key,
                               "G")
indicsp$run_indicator_analysis(dataset_name,
                               tables[["raw_clade"]],
                               ind_sp_dir_clade,
                               treatment_key,
                               "G")

# Relative Abundance ------------------------------------------------------
exploratory$make_interest_abundance(tables[["raw_clade"]],
                                    treatment_key,
                                    dataset_name,
                                    additional_taxa,
                                    rel_abund_dir)


# Alpha Diversity ---------------------------------------------------------
exploratory$make_all_alpha_plots(tables[["raw_clade"]],
                                 treatment_key,
                                 dataset_name,
                                 alpha_div_dir)


# Beta Diversity ----------------------------------------------------------
exploratory$make_nmds_plots(tables[["raw_clade"]],
                            treatment_key,
                            dataset_name,
                            nmds_dir)
# _________________________________________________________________________
