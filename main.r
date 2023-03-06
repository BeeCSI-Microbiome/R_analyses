# Author(s): Lance Lansing, Jonathan Ho, Kurt Clarke

# Main script for analysis of taxonomic classification results


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


run_all_analyses <- function(dataset_name, counts_path, treatment_key,
                             filter_string="", act1_meta=NULL, outdir="results/"){
  # Globals -----------------------------------------------------------------
  # Name of the dataset for file writing purposes
  # dataset_name <- "a_flx_d1_low_exposure"
  # filepath to the taxon count table
  # counts_path <- "../data/a_flx_2021/a_flx_2021_aggregated_counts.csv"

  # Activity 1 Metadata sheet
  # act1_meta <- read_csv("activity-1_metadata.csv")


  # Read input files --------------------------------------------------------
  ct <- read_csv(counts_path)

  # filter the treatment that isn't being analyzed
  if (filter_string != "")
    ct <- ct %>% select(-matches(filter_string))
  # _________________________________________________________________________

  # The following key should match the substring in the sample name that specifies
  # a treatment with the name of that treatment, including control
  # THE CONTROL TREATMENT MUST BE THE FIRST ELEMENT (e.g. control, unexposed)
  # treatment_key <- list(d0="control",d1="low_exposure")

  ## Create output directories ####
  main_outdir <- glue("{outdir}/{dataset_name}")
  nmds_dir <- glue("{main_outdir}/nmds_anosim")
  alpha_div_dir <- glue("{main_outdir}/alpha_diversity")
  rel_abund_dir <-  glue("{main_outdir}/relative_abundance")
  da_dir <- glue("{main_outdir}/differential_abundance")
  da_ancombc_dir <- glue("{da_dir}/ancombc")
  ind_sp_dir <- glue("{da_dir}/indicator_species_analysis")

  create_dir_if_nonexistant <- function(path) {
    ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
  }

  lapply(c(outdir, main_outdir, nmds_dir, alpha_div_dir, rel_abund_dir, da_dir,
           da_ancombc_dir,ind_sp_dir), create_dir_if_nonexistant)
  # _________________________________________________________________________


  # Reads summary -----------------------------------------------------------
  # Creates a table containing total, unclassified, and classified read counts,
  # and % of reads classified as A. mellifera, Bacteria, and other specific taxa
  reads_summary <- rsummary$create_summary_table(ct, dataset_name, main_outdir)
  # _________________________________________________________________________

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
  ancom_species_results <- da_ancombc$run_ancombc(tables[["raw_clade"]],
                                                  treatment_key,
                                                  dataset_name,
                                                  "S",
                                                  da_ancombc_dir,
                                                  act1_meta)
  ancom_genus_results <- da_ancombc$run_ancombc(tables[["raw_clade"]],
                                                treatment_key,
                                                dataset_name,
                                                "G",
                                                da_ancombc_dir,
                                                act1_meta)

  # Indicator Taxa Analysis--------------------------------------------------
  # indicsp$run_indicator_analysis(tables[["raw_clade"]],
  #                                treatment_key,
  #                                dataset_name,
  #                                ind_sp_dir)
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

  # Other analyses ----------------------------------------------------------

  # Get additional taxa above 1% rel. abundance
  taxa_cutoff_table <- exploratory$taxa_cutoff_explore(tables[["raw_clade"]], dataset_name)
  additional_taxa <- taxa_cutoff_table$name[taxa_cutoff_table$new_taxa_of_interest]

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

  # Write the pdf report ----------------------------------------------------
  rmarkdown::render("main.qmd",
                    output_file = glue("{main_outdir}/{dataset_name}_report_{Sys.Date()}.pdf"),
                    params = list(dataset_name = dataset_name))
}
