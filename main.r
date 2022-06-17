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
widen_results <- use("scripts/widenResults.R")
da_fitzig <- use("scripts/da_fitzig.R")
da_ancombc <- use("scripts/da-ancombc.R")
# _________________________________________________________________________


# Globals -----------------------------------------------------------------
# Name of the dataset for file writing purposes
dataset_name <- "oxy_2021"

# Create strings for output directories
main_outdir <- glue("results/{dataset_name}")
nmds_dir <- glue("{main_outdir}/nmds_anosim")
alpha_div_dir <- glue("{main_outdir}/alpha_diversity")
rel_abund_dir <-  glue("{main_outdir}/relative_abundance")
# <DA> Subdirectory for kraken report matrices
kraken_matrix_dir <-
  glue("{main_outdir}/aggregated_kraken_reports")
# <DA> Subdirectory for DA results
da_dir <- glue("{main_outdir}/differential_abundance")
da_fitzig_dir <- glue("{da_dir}/fitzig")
da_ancombc_dir <- glue("{da_dir}/ancombc")

create_dir_if_nonexistant <- function(path) {
  ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
}
## Create output directories ####
lapply(c(main_outdir, nmds_dir, alpha_div_dir, rel_abund_dir,
         kraken_matrix_dir, da_dir, da_fitzig_dir, da_ancombc_dir),
       create_dir_if_nonexistant)

## Input file paths ####
# Counts (clade and taxon counts, uncollapsed) and metadata tables
counts_path <- "../data/oxy_2021/oxy_2021_all_taxa.csv"
metadata_filepath <- "../data/oxy_2021/oxy_2021_metadata.csv"

# Path to directory with kraken reports
kraken_report_dir <- "../data/oxy_2021/kraken_reports"
# Kraken report paths for fitzig DA
krakenReportPaths <-
  Sys.glob(glue("{kraken_report_dir}/*_report.txt"))
  # c(Sys.glob(glue("{kraken_report_dir}/*dT*_report.txt")))

krakenReportNames <- list.files(path = kraken_report_dir,
                                pattern = "*_report.txt")

## User-defined values ####
# Dataset-specific taxa of interest (provide a list of taxa strings)
# e.g. c("Lactobacillus", "Gilliamella apis")
additional_taxa <- c("Pantoea agglomerans")


# The following key should match the substring in the sample name that specifies 
# a treatment with the name of that treatment, including control
treatment_key <- list(d0="control",d1="oxytetracycline")
# a treatment with the name of that treatment
# THE CONTROL TREATMENT MUST BE THE FIRST ELEMENT (e.g. control, unexposed)

# Percentile value used by CSS (default = 0.5)
css_percentile <- 0.5

# # Statistical analyses list for DA
# statistical_analyses <- list(
#   # ACTIVITY 2 - Example: Corn (likely applicable to most activity 2 datasets)
#   # Description: Compare exposed vs. unexposed, use replicate as random effect
#   list(
#     name = "Exposure",
#     subsets = list(),
#     model_matrix = "~ 0 + Exposed",
#     contrasts = list(
#       "ExposedTRUE - ExposedFALSE"
#     ),
#     random_effect = "Replicate"
#   )
# )

statistical_analyses = list(
  # ACTIVITY 1 - Example: CTX
  # Description: Compare control and 2 treatments, use replicate as random effect
  list(
    name = "Treatment",
    subsets = list(),
    model_matrix = "~ 0 + Treatment",
    contrasts = list(
      "Treatmentcontrol - Treatmentoxytetracycline"
    ),
    random_effect = "Replicate"
  )
)
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


# Differential Abundance --------------------------------------------------
widen_results$widen_results_function(krakenReportPaths,
                                     krakenReportNames,
                                     kraken_matrix_dir)
da_fitzig$kraken_differential_abundance(
  kraken_matrix_dir,
  metadata_filepath,
  da_fitzig_dir,
  statistical_analyses
)


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
