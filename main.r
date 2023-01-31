# Author(s): Lance Lansing, Jonathan Ho, Kurt Clarke

# Main script for analysis of taxonomic classification results

# Input: table from Pavian with clades AND taxon counts, not collapsed
# User defined inputs (eg input file paths) must be set under section "Globals"

# setwd("C:/Users/clarkeku/OneDrive - AGR-AGR/Documents/GitHub/R_analyses-main (4)/R_analyses-main")

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


run_all_analyses <- function(dataset_name, counts_path, treatment_key,  filter_string=""){
  # Globals -----------------------------------------------------------------
  # Name of the dataset for file writing purposes
  # dataset_name <- "a_flx_d1_low_exposure"
  # filepath to the taxon count table
  # counts_path <- "../data/a_flx_2021/a_flx_2021_aggregated_counts.csv"

  # Activity 1 Metadata sheet
  # act1_meta <- read_csv("activity-1_metadata.csv")
  act1_meta <- NULL

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
  main_outdir <- glue("results/act2_results/{dataset_name}")
  nmds_dir <- glue("{main_outdir}/nmds_anosim")
  alpha_div_dir <- glue("{main_outdir}/alpha_diversity")
  rel_abund_dir <-  glue("{main_outdir}/relative_abundance")
  da_dir <- glue("{main_outdir}/differential_abundance")
  da_ancombc_dir <- glue("{da_dir}/ancombc")
  ind_sp_dir <- glue("{da_dir}/indicator_species_analysis")

  create_dir_if_nonexistant <- function(path) {
    ifelse(!dir.exists(path), dir.create(path, mode = "777"), FALSE)
  }

  lapply(c(main_outdir, nmds_dir, alpha_div_dir, rel_abund_dir, da_dir,
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

# Activity 2
# SOY_2020
run_all_analyses(dataset_name = "soy_2020",
                 counts_path = "../data/soy_2020/soy_2020_aggregated_counts.csv",
                 treatment_key = list(u="unexposed",e="exposed"))


 # Activity 1
# # A_BOS
# run_all_analyses(dataset_name = "a_bos_d1_low_exposure",
#                  counts_path = "../data/a_bos_2021/a_bos_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "a_bos_d2_high_exposure",
#                  counts_path = "../data/a_bos_2021/a_bos_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
# # A_FLX
# run_all_analyses(dataset_name = "a_flx_d1_low_exposure",
#                  counts_path = "../data/a_flx_2021/a_flx_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "a_flx_d2_high_exposure",
#                  counts_path = "../data/a_flx_2021/a_flx_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
# # A_PYM
# run_all_analyses(dataset_name = "a_pym_d1_low_exposure",
#                  counts_path = "../data/a_pym_2021/a_pym_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "a_pym_d2_high_exposure",
#                  counts_path = "../data/a_pym_2021/a_pym_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
#
#
# # B_FLY
# run_all_analyses(dataset_name = "b_fly_d1_low_exposure",
#                  counts_path = "../data/b_fly_2021/b_fly_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "b_fly_d2_high_exposure",
#                  counts_path = "../data/b_fly_2021/b_fly_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
# # B_PYC
# run_all_analyses(dataset_name = "b_pyc_d1_low_exposure",
#                  counts_path = "../data/b_pyc_2021/b_pyc_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "b_pyc_d2_high_exposure",
#                  counts_path = "../data/b_pyc_2021/b_pyc_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
#
# # c_chl
# run_all_analyses(dataset_name = "c_chl_d1_low_exposure",
#                  counts_path = "../data/c_chl_2021/c_chl_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "c_chl_d2_high_exposure",
#                  counts_path = "../data/c_chl_2021/c_chl_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
# # c_spn
# run_all_analyses(dataset_name = "c_spn_d1_low_exposure",
#                  counts_path = "../data/c_spn_2021/c_spn_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "c_spn_d2_high_exposure",
#                  counts_path = "../data/c_spn_2021/c_spn_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
# # c_spr
# run_all_analyses(dataset_name = "c_spr_d1_low_exposure",
#                  counts_path = "../data/c_spr_2021/c_spr_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "c_spr_d2_high_exposure",
#                  counts_path = "../data/c_spr_2021/c_spr_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
#
# # d_flp
# run_all_analyses(dataset_name = "d_flp_d1_low_exposure",
#                  counts_path = "../data/d_flp_2021/d_flp_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "d_flp_d2_high_exposure",
#                  counts_path = "../data/d_flp_2021/d_flp_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
# # d_sul
# run_all_analyses(dataset_name = "d_sul_d1_low_exposure",
#                  counts_path = "../data/d_sul_2021/d_sul_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "d_sul_d2_high_exposure",
#                  counts_path = "../data/d_sul_2021/d_sul_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
#
# # e_gly
# run_all_analyses(dataset_name = "e_gly_d1_low_exposure",
#                  counts_path = "../data/e_gly_2021/e_gly_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "e_gly_d2_high_exposure",
#                  counts_path = "../data/e_gly_2021/e_gly_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
# # e_met
# run_all_analyses(dataset_name = "e_met_d1_low_exposure",
#                  counts_path = "../data/e_met_2021/e_met_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_exposure"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "e_met_d2_high_exposure",
#                  counts_path = "../data/e_met_2021/e_met_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_exposure"),
#                  filter_string = "_d1")
#
# # CLO
# run_all_analyses(dataset_name = "clo_2020_d1_sublethal",
#                  counts_path = "../data/clo_2020/clo_2020_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="sublethal"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "clo_2020_d2_acute",
#                  counts_path = "../data/clo_2020/clo_2020_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="acute"),
#                  filter_string = "_d1")
# # THI
# run_all_analyses(dataset_name = "thi_2020_d1_sublethal",
#                  counts_path = "../data/thi_2020/thi_2020_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="sublethal"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "thi_2020_d2_acute",
#                  counts_path = "../data/thi_2020/thi_2020_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="acute"),
#                  filter_string = "_d1")
# # CTX
# run_all_analyses(dataset_name = "ctx_2020_dC_clothianidin",
#                  counts_path = "../data/ctx_2020/ctx_dC_2020_aggregated_counts.csv",
#                  treatment_key = list(d0="control",dC="sublethal"))
# run_all_analyses(dataset_name = "ctx_2020_dT_thiamethoxam",
#                  counts_path = "../data/ctx_2020/ctx_dT_2020_aggregated_counts.csv",
#                  treatment_key = list(d0="control",dT="acute"))
# # CFS
# run_all_analyses(dataset_name = "cfs_2021_dC_chlorantraniliprole",
#                  counts_path = "../data/cfs_dC_2021/cfs_dC_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",dC="chlorantraniliprole"))
# run_all_analyses(dataset_name = "cfs_2021_dF_flupyradifurone",
#                  counts_path = "../data/cfs_dF_2021/cfs_dF_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",dF="flupyradifurone"))
# run_all_analyses(dataset_name = "cfs_2021_dS_sulfoxaflor",
#                  counts_path = "../data/cfs_dS_2021/cfs_dS_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",dS="sulfoxaflor"))
#
# # AFB
# run_all_analyses(dataset_name = "afb_2021_d1_subclinical",
#                  counts_path = "../data/afb_2021/afb_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="subclinical"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "afb_2021_d2_clinical",
#                  counts_path = "../data/afb_2021/afb_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="clinical"),
#                  filter_string = "_d1")
#
# # IAP
# run_all_analyses(dataset_name = "iap_2021",
#                  counts_path = "../data/iap_2021/iap_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="infected"))
#
# # NOS
# run_all_analyses(dataset_name = "nos_2021_d1_low_dose",
#                  counts_path = "../data/nos_2021/nos_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_dose"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "nos_2021_d2_high_dose",
#                  counts_path = "../data/nos_2021/nos_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_dose"),
#                  filter_string = "_d1")
#
# # VAR
# run_all_analyses(dataset_name = "var_2021_d1_low_dose",
#                  counts_path = "../data/var_2021/var_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_dose"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "var_2021_d2_high_dose",
#                  counts_path = "../data/var_2021/var_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_dose"),
#                  filter_string = "_d1")
# # OXY
# run_all_analyses(dataset_name = "oxy_2021",
#                  counts_path = "../data/oxy_2021/oxy_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="oxytetracycline"))
#
# # PDV
# run_all_analyses(dataset_name = "pdv_2021_d1_nutritious_monofloral",
#                  counts_path = "../data/pdv_2021/pdv_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="nutritious_monofloral"),
#                  filter_string = "_d2|_t1")
# run_all_analyses(dataset_name = "pdv_2021_d2_non_nutritious_polyfloral",
#                  counts_path = "../data/pdv_2021/pdv_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="non_nutritious_polyfloral"),
#                  filter_string = "_d1|_t1")
#
# # PRE
# run_all_analyses(dataset_name = "pre_2021_d1_restricted",
#                  counts_path = "../data/pre_2021/pre_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="restricted"),
#                  filter_string = "_d2|_t1")
# run_all_analyses(dataset_name = "pre_2021_d2_supplemented",
#                  counts_path = "../data/pre_2021/pre_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="supplemented"),
#                  filter_string = "_d1|_t1")
#
# # AMZ
# run_all_analyses(dataset_name = "amz_2021",
#                  counts_path = "../data/amz_2021/amz_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="amitraz"))
# # OXA
# run_all_analyses(dataset_name = "oxa_2021",
#                  counts_path = "../data/oxa_2021/oxa_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="oxalic_acid"))
# # CHA
# run_all_analyses(dataset_name = "cha_2021",
#                  counts_path = "../data/cha_2021/cha_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="chalkbrood"))
#
# # FMA
# run_all_analyses(dataset_name = "fma_2021_d1_low_dose",
#                  counts_path = "../data/fma_2021/fma_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d1="low_dose"),
#                  filter_string = "_d2")
# run_all_analyses(dataset_name = "fma_2021_d2_high_dose",
#                  counts_path = "../data/fma_2021/fma_2021_aggregated_counts.csv",
#                  treatment_key = list(d0="control",d2="high_dose"),
#                  filter_string = "_d1")
#
