# Author(s): Lance Lansing, Jonathan Ho, Kurt Clarke

# Main script for analysis of taxonomic classification results

# Input: table from Pavian with clades AND taxon counts, not collapsed
# User defined inputs (eg input file paths) must be set under section "Globals"

setwd("C:/Users/clarkeku/OneDrive - AGR-AGR/Documents/GitHub/R_analyses-main (4)/R_analyses-main")

# Package setup -----------------------------------------------------------
packages <- c("tidyverse",
              "vegan",
              "modules",
              "data.table",
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
# filepath to the taxon count table
counts_path <- "C:/Users/clarkeku/OneDrive - AGR-AGR/Documents/GitHub/R_analyses-main (4)/R_analyses-main/data/oxy_2021/oxy_2021_aggregated_counts.csv"

# Create strings for output directories
main_outdir <- glue("results/{dataset_name}")
nmds_dir <- glue("{main_outdir}/nmds_anosim")
alpha_div_dir <- glue("{main_outdir}/alpha_diversity")
rel_abund_dir <-  glue("{main_outdir}/relative_abundance")
da_dir <- glue("{main_outdir}/differential_abundance")
da_ancombc_dir <- glue("{da_dir}/ancombc")
ind_sp_dir <- glue("{da_dir}/indicator_species_analysis")

create_dir_if_nonexistant <- function(path) {
  ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
}
## Create output directories ####
lapply(c(main_outdir, nmds_dir, alpha_div_dir, rel_abund_dir, da_dir,
         da_ancombc_dir,ind_sp_dir), create_dir_if_nonexistant)


## User-defined values ####
# Dataset-specific taxa of interest (provide a list of taxa strings)
# e.g. c("Lactobacillus", "Gilliamella apis")
additional_taxa <- c(cor_2020=NA,
                     cac_2020=NA,
                     cas_2020=NA,
                     cra_2020=NA,
                     hbb_2020_t1=NA,
                     hbb_2020_t2=c("Spiroplasma melliferum",
                                   "Paenibacillus alvei"),
                     hbb_2020_t3=c("Spiroplasma melliferum",
                                   "Paenibacillus alvei"),
                     hbb_2020_t4=c("Paenibacillus alvei"),
                     soy_2020=c("Pantoea agglomerans"),
                     clo_2020=c("Spiroplasma melliferum",
                                "Serratia marcescens"),
                     ctx_2020_dC=c("Pantoea agglomerans"),
                     ctx_2020_dT=c("Pantoea agglomerans"),
                     thi_2020=c("Spiroplasma apis"),
                     a_bos_2021=NA,
                     a_flx_2021=NA,
                     a_pym_2021=NA,
                     app_2021=NA,
                     b_fly_2021=NA,
                     b_pyc_2021=NA,
                     c_chl_2021=NA,
                     c_spn_2021=NA,
                     c_spr_2021=NA,
                     cac_2021=NA,
                     cas_2021=c("Spiroplasma melliferum"),
                     cfs_dC_2021=c("Pantoea agglomerans"),
                     cfs_dF_2021=c("Pantoea agglomerans"),
                     cfs_dS_2021=c("Pantoea agglomerans"),
                     cra_2021=c("Paenibacillus alvei"),
                     d_flp_2021=NA,
                     d_sul_2021=NA,
                     e_gly_2021=c("Serratia marcescens"),
                     e_met_2021=NA,
                     hbb_2021_t1=NA,
                     hbb_2021_t2=c("Spiroplasma melliferum",
                                   "Paenibacillus alvei"),
                     hbb_2021_t3=c("Spiroplasma melliferum",
                                   "Paenibacillus alvei"),
                     hbb_2021_t4=c("Spiroplasma melliferum",
                                   "Paenibacillus alvei"),
                     lbb_2021=NA,
                     pdv_2021_t1=NA,
                     pdv_2021_t2=NA,
                     pre_2021_t1=NA,
                     pre_2021_t2=NA,
                     afb_2021=NA,
                     cha_2021=NA,
                     iap_2021=c("Bombella intestini"),
                     nos_2021=NA,
                     var_2021=c("Bombella intestini"))
additional_taxa <- additional_taxa[dataset_name]


# The following key should match the substring in the sample name that specifies 
# a treatment with the name of that treatment, including control
treatment_key <- list(d0="control",d1="oxytetracycline")

# treatment key for generalized .qmd document output. 
treatmentkey.df = as.data.frame(do.call(cbind, treatment_key))
write.csv(treatmentkey.df, glue("{main_outdir}/treatmentKey.csv"),
          row.names = FALSE)

# a treatment with the name of that treatment
# THE CONTROL TREATMENT MUST BE THE FIRST ELEMENT (e.g. control, unexposed)
# _________________________________________________________________________


# Read input files --------------------------------------------------------
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
indicsp$run_indicator_analysis(tables[["raw_clade"]],
                               treatment_key,
                               dataset_name,
                               ind_sp_dir)
indicsp$run_indicator_analysis(tables[["raw_clade"]],
                               treatment_key,
                               dataset_name,
                               ind_sp_dir,
                               "S")
indicsp$run_indicator_analysis(tables[["raw_clade"]],
                               treatment_key,
                               dataset_name,
                               ind_sp_dir,
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

