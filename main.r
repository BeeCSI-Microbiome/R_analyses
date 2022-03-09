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
scaling <- use("scripts/reads_scaling.R")
exploratory <- use("scripts/exploratory_functions.R")
widen_results <- use("scripts/widenResults.R")
da <- use("scripts/differential_abundance.R")
# _________________________________________________________________________


# Globals -----------------------------------------------------------------
# Name of the dataset for file writing purposes
dataset_name <- "ctx_2020"

# Create strings for output directories
main_outdir <- glue("results/{dataset_name}")
nmds_dir <- glue("{main_outdir}/nmds_anosim")
alpha_div_dir <- glue("{main_outdir}/alpha_diversity")
rel_abund_dir <-  glue("{main_outdir}/relative_abundance")
# <DA> Subdirectory for kraken report matrices
kraken_matrix_dir <- glue("{main_outdir}/aggregated_kraken_reports")
# <DA> Subdirectory for DA results
da_dir <- glue("{main_outdir}/differential_abundance")

## Create output directories ####
ifelse(!dir.exists(main_outdir), dir.create(main_outdir, mode = "777"), FALSE)
ifelse(!dir.exists(nmds_dir), dir.create(nmds_dir, mode = "777"), FALSE)
ifelse(!dir.exists(alpha_div_dir), dir.create(alpha_div_dir, mode = "777"), FALSE)
ifelse(!dir.exists(rel_abund_dir), dir.create(rel_abund_dir, mode = "777"), FALSE)
ifelse(!dir.exists(kraken_matrix_dir), dir.create(kraken_matrix_dir, mode = "777"), FALSE)
ifelse(!dir.exists(da_dir), dir.create((da_dir), mode = "777"), FALSE)



## Input file paths ####
# Counts (clade and taxon counts, uncollapsed) and metadata tables
counts_path <- "../data/cas_2020/cas_all_taxa.tsv"
metadata_filepath <- "../data/cas_2020/cas_2020_metadata.csv"

# Path to directory with kraken reports
kraken_report_dir <- "../data/ctx_2020/kraken_reports"
# Kraken report paths <DA>
krakenReportPaths <- Sys.glob(glue("{kraken_report_dir}/*_report.txt"))
krakenReportNames <- list.files(path = kraken_report_dir,
                                pattern = "*_report.txt")

## User-defined values ####
# Dataset-specific taxa of interest (provide a list of taxa strings)
# e.g. c("Lactobacillus", "Gilliamella apis")
additional_taxa <- NA

# Treatment and replicate names
treat_names <- c("control", "clothianidin", "thiamethoxam")
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")

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
      "Treatmentcontrol - Treatmentclothianidin",
      "Treatmentcontrol - Treatmentthiamethoxam",
      "Treatmentclothianidin - Treatmentthiamethoxam"
    ),
    random_effect = "Replicate"
  )
)
# _________________________________________________________________________


# Read input files --------------------------------------------------------
# Count table
ct <- read_tsv(counts_path)
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
da$kraken_differential_abundance(kraken_matrix_dir,
                                 metadata_filepath,
                                 da_dir,
                                 statistical_analyses,
                                 css_percentile)


# Formatting, filtering, and scaling --------------------------------------
# Format and perform filtering on the table
ct <- ip$format_count_table(ct) %>%
  ip$filter_table() %>%
  ip$group_taxa_of_interest()

# Return list of MRexperiment object, and raw and scaled taxon and clade tables
tables <- scaling$scaling_procedure(ct, css_percentile, dataset_name)
# __________________________________________________________________________


# Relative Abundance ------------------------------------------------------
exploratory$make_interest_abundance(tables[["raw_clade"]],
                                    treat_names,
                                    rep_names,
                                    dataset_name,
                                    additional_taxa,
                                    rel_abund_dir)


# Alpha Diversity ---------------------------------------------------------
exploratory$make_all_alpha_plots(tables[["raw_clade"]],
                                 treat_names,
                                 rep_names,
                                 dataset_name,
                                 alpha_div_dir)


# Beta Diversity ----------------------------------------------------------
exploratory$make_nmds_plots(tables[["raw_clade"]],
                            treat_names,
                            rep_names,
                            dataset_name,
                            nmds_dir)
# _________________________________________________________________________

