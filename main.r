# Author(s): Lance Lansing, Jonathan Ho, Kurt Clarke

# Main script for analysis of taxonomic classification results

# Input: table from Pavian with clades AND taxon counts, not collapsed
# User defined inputs (eg input file paths) must be set under section "Globals"

setwd("C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses")

# Package setup -----------------------------------------------------------
packages <- c("tidyverse",
              "vegan",
              "modules",
              "data.table",
              "ggplot2",
              "glue")
lapply(packages, library, character.only = TRUE)
library(dplyr)
# _________________________________________________________________________


# Load aux scripts as modules ---------------------------------------------
# rsummary <- use("scripts/reads_summary.R")
# ip <- use("scripts/initial_processing.R")
exploratory <- use("scripts/exploratory_functions.R")
# da_ancombc <- use("scripts/da-ancombc.R")
# indicsp <- use("scripts/indicator_taxa_analysis.R")
# _________________________________________________________________________


# Globals -----------------------------------------------------------------
# Name of the dataset for file writing purposes
dataset_name <- "1beta_acrosscrops_allprov_withhbb21_onlyexposed"
# filepath to the taxon count table
counts_path <- "C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_app21.csv"

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
# additional_taxa <- c(cor_2020=NA,
#                      cac_2020=NA,
#                      cas_2020=NA,
#                      cra_2020=NA,
#                      hbb_2020_t1=NA,
#                      hbb_2020_t2=c("Spiroplasma melliferum",
#                                    "Paenibacillus alvei"),
#                      hbb_2020_t3=c("Spiroplasma melliferum",
#                                    "Paenibacillus alvei"),
#                      hbb_2020_t4=c("Paenibacillus alvei"),
#                      soy_2020=c("Pantoea agglomerans"),
#                      clo_2020=c("Spiroplasma melliferum",
#                                 "Serratia marcescens"),
#                      ctx_2020_dC=c("Pantoea agglomerans"),
#                      ctx_2020_dT=c("Pantoea agglomerans"),
#                      thi_2020=c("Spiroplasma apis"),
#                      a_bos_2021=NA,
#                      a_flx_2021=NA,
#                      a_pym_2021=NA,
#                      app_2021=NA,
#                      b_fly_2021=NA,
#                      b_pyc_2021=NA,
#                      c_chl_2021=NA,
#                      c_spn_2021=NA,
#                      c_spr_2021=NA,
#                      cac_2021=NA,
#                      cas_2021=c("Spiroplasma melliferum"),
#                      cfs_dC_2021=c("Pantoea agglomerans"),
#                      cfs_dF_2021=c("Pantoea agglomerans"),
#                      cfs_dS_2021=c("Pantoea agglomerans"),
#                      cra_2021=c("Paenibacillus alvei"),
#                      d_flp_2021=NA,
#                      d_sul_2021=NA,
#                      e_gly_2021=c("Serratia marcescens"),
#                      e_met_2021=NA,
#                      hbb_2021_t1=NA,
#                      hbb_2021_t2=c("Spiroplasma melliferum",
#                                    "Paenibacillus alvei"),
#                      hbb_2021_t3=c("Spiroplasma melliferum",
#                                    "Paenibacillus alvei"),
#                      hbb_2021_t4=c("Spiroplasma melliferum",
#                                    "Paenibacillus alvei"),
#                      lbb_2021=NA,
#                      pdv_2021_t1=NA,
#                      pdv_2021_t2=NA,
#                      pre_2021_t1=NA,
#                      pre_2021_t2=NA,
#                      afb_2021=NA,
#                      cha_2021=NA,
#                      iap_2021=c("Bombella intestini"),
#                      nos_2021=NA,
#                      var_2021=c("Bombella intestini"))
# additional_taxa <- additional_taxa[dataset_name]


# The following key should match the substring in the sample name that specifies 
# a treatment with the name of that treatment, including control
treatment_key <- c(u="unexposed",e="exposed")
# a treatment with the name of that treatment
# THE CONTROL TREATMENT MUST BE THE FIRST ELEMENT (e.g. control, unexposed)
# _________________________________________________________________________


# Read and format primary crop input file --------------------------------------------------------
ct <- read_csv(counts_path)

ct = ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(ct) = paste(colnames(ct), str_sub(counts_path, -6, -5), sep="_")
ct <- select(ct, !(matches("t1|t3|t4")))
ct <- separate(ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Read and format all other crop input files________________________________________________________
#Crop 2-----
crop2_counts_path = "C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_cac20.csv"

crop2_ct = read_csv(crop2_counts_path)
crop2_ct = crop2_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop2_ct) = paste(colnames(crop2_ct), str_sub(crop2_counts_path, -6, -5), sep="_")
crop2_ct <- select(crop2_ct, !(matches("t1|t3|t4")))
crop2_ct <- separate(crop2_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 3-----
crop3_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_cas20.csv"

crop3_ct = read_csv(crop3_counts_path)
crop3_ct = crop3_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop3_ct) = paste(colnames(crop3_ct), str_sub(crop3_counts_path, -6, -5), sep="_")
crop3_ct <- select(crop3_ct, !(matches("t1|t3|t4")))
crop3_ct <- separate(crop3_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 4------
crop4_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_cac21.csv"

crop4_ct = read_csv(crop4_counts_path)
crop4_ct = crop4_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop4_ct) = paste(colnames(crop4_ct), str_sub(crop4_counts_path, -6, -5), sep="_")
crop4_ct <- select(crop4_ct, !(matches("t1|t3|t4")))
crop4_ct <- separate(crop4_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 5------
crop5_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_cas21.csv"

crop5_ct = read_csv(crop5_counts_path)
crop5_ct = crop5_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop5_ct) = paste(colnames(crop5_ct), str_sub(crop5_counts_path, -6, -5), sep="_")
crop5_ct <- select(crop5_ct, !(matches("t1|t3|t4")))
crop5_ct <- separate(crop5_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 6------
crop6_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_cor20.csv"

crop6_ct = read_csv(crop6_counts_path)
crop6_ct = crop6_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop6_ct) = paste(colnames(crop6_ct), str_sub(crop6_counts_path, -6, -5), sep="_")
crop6_ct <- select(crop6_ct, !(matches("t1|t3|t4")))
crop6_ct <- separate(crop6_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 7------
crop7_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_cra20.csv"

crop7_ct = read_csv(crop7_counts_path)
crop7_ct = crop7_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop7_ct) = paste(colnames(crop7_ct), str_sub(crop7_counts_path, -6, -5), sep="_")
crop7_ct <- select(crop7_ct, !(matches("t1|t3|t4")))
crop7_ct <- separate(crop7_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 8------
crop8_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_cra21.csv"

crop8_ct = read_csv(crop8_counts_path)
crop8_ct = crop8_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop8_ct) = paste(colnames(crop8_ct), str_sub(crop8_counts_path, -6, -5), sep="_")
crop8_ct <- select(crop8_ct, !(matches("t1|t3|t4")))
crop8_ct <- separate(crop8_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 9------
crop9_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_hbb_t2_20.csv"

crop9_ct = read_csv(crop9_counts_path)
crop9_ct = crop9_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop9_ct) = paste(colnames(crop9_ct), str_sub(crop9_counts_path, -6, -5), sep="_")
crop9_ct <- select(crop9_ct, !(matches("t1|t3|t4")))
crop9_ct <- separate(crop9_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 10------
crop10_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_hbb21.csv"

crop10_ct = read_csv(crop10_counts_path)
crop10_ct = crop10_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop10_ct) = paste(colnames(crop10_ct), str_sub(crop10_counts_path, -6, -5), sep="_")
crop10_ct <- select(crop10_ct, !(matches("t1|t3|t4")))
crop10_ct <- separate(crop10_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 11------
crop11_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_lbb21.csv"

crop11_ct = read_csv(crop11_counts_path)
crop11_ct = crop11_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop11_ct) = paste(colnames(crop11_ct), str_sub(crop11_counts_path, -6, -5), sep="_")
crop11_ct <- select(crop11_ct, !(matches("t1|t3|t4")))
crop11_ct <- separate(crop11_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#Crop 12------
crop12_counts_path ="C:/Users/obrienj/OneDrive - AGR-AGR/Desktop/Bee Stuff/R_analyses/AMR_crop_data/AMR_analytic_matrix_soy20.csv"

crop12_ct = read_csv(crop12_counts_path)
crop12_ct = crop12_ct %>% 
  filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
colnames(crop12_ct) = paste(colnames(crop12_ct), str_sub(crop12_counts_path, -6, -5), sep="_")
crop12_ct <- select(crop12_ct, !(matches("t1|t3|t4")))
crop12_ct <- separate(crop12_ct, c(1), into = c("geneID", "type", "class", "mechanism", "gene"), sep = "/")

#'*reads summary actually looks useful. Lance recommended to skip for now*
# Reads summary -----------------------------------------------------------
# Creates a table containing total, unclassified, and classified read counts,
# and % of reads classified as A. mellifera, Bacteria, and other specific taxa
#reads_summary <- rsummary$create_summary_table(ct, dataset_name)
# _________________________________________________________________________

merged_crops = merge(ct, crop2_ct, all = TRUE) %>% 
  merge(crop3_ct, all = TRUE) %>% 
  merge(crop4_ct, all = TRUE)%>% 
  merge(crop5_ct, all = TRUE) %>% 
  merge(crop6_ct, all = TRUE)%>% 
  merge(crop7_ct, all = TRUE) %>% 
  merge(crop8_ct, all = TRUE) %>% 
  merge(crop9_ct, all = TRUE)%>% 
  merge(crop10_ct, all = TRUE) %>%
  merge(crop11_ct, all = TRUE)%>% 
  merge(crop12_ct, all = TRUE)
  
# merged_crops = merged_crops %>% 
#   filter(!str_detect(gene_accession, "RequiresSNPConfirmation"))
# 
# merged_crops <- select(merged_crops, !(matches("t1|t3|t4")))
# merged_crops <- separate(merged_crops, col = gene_accession, into = c("geneID", "type", "class", "mechanism", "gene", "SNP"), sep = "/")

# Beta Diversity All Crops----------------------------------------------------------
#'*replaced by across crops.  other functions call some of the "all" functions, so can't delete them*
# exploratory$make_nmds_plots_all_crops(merged_crops,
#                             treatment_key,
#                             dataset_name,
#                             nmds_dir)

#'*requested in meeting*
vIew_plot = exploratory$make_nmds_plots_across_crops(merged_crops,
                                      treatment_key,
                                      dataset_name,
                                      nmds_dir)
# _________________________________________________________________________

# Beta Diversity Between All Provinces----------------------------------------------------------
ploto = exploratory$make_nmds_plots_provinces(merged_crops,
                                      treatment_key,
                                      dataset_name,
                                      nmds_dir)
# _________________________________________________________________________
#'*requested in meeting*
# Beta Diversity for All Crops of a Province, Between Paired Provinces.  Unexposed crops lumped together as one "crop"---
paired = exploratory$make_nmds_paired_prov_unex_by_crop(merged_crops,
                                                      treatment_key,
                                                      dataset_name,
                                                      nmds_dir)

# Beta Diversity comparing Crops between years.  Each crop looked at individually.  Unexposed and exposed together
exploratory$make_nmds_paired_crop_by_year(merged_crops,
                                                        treatment_key,
                                                        dataset_name,
                                                        nmds_dir)


# Beta Diversity Between Paired Provinces (WIP)----------------------------------------------------------
exploratory$make_nmds_plots_paired_provinces(merged_crops,
                                              treatment_key,
                                              dataset_name,
                                              nmds_dir)
# _________________________________________________________________________

# Differential Abundance --------------------------------------------------

# Formatting, filtering, and calculating clade counts ---------------------
 #Format and perform filtering on the table
#'* How much of ip do we actually need here?  Don't need to get rid of taxLineage, clade reads, etc. Is it changing the format from csv?*
#''*do need "tables" - that's what gets fed to the RA scripts*
#''*do NOT need format_count_table - the formatting I wrote above (ct <- separate)) is doing that*
 # ct <- ip$format_count_table(ct) %>%
 #   ip$filter_table() %>%
 #   ip$group_taxa_of_interest()

#'*I think tables is more or less the same format as ct already is, so perhaps can just sub in ct instead of tables?*

#'  tables <- ip$calculate_clade_counts(ct)
#' 
#' #'*putting the newly separated clade and taxon counts into their own separate csv's*
#'  write.csv(tables[["raw_clade"]],
#'            glue("{main_outdir}/raw_clade_counts.csv"),
#'            row.names = FALSE)
#'  write.csv(tables[["raw_taxon"]],
#'            glue("{main_outdir}/raw_taxon_counts.csv"),
#'            row.names = FALSE)


# # ANCOMBC Differential Abundance ------------------------------------------
# da_ancombc$run_ancombc(tables[["raw_clade"]],
#                        treatment_key,
#                        dataset_name,
#                        "S",
#                        da_ancombc_dir)
# da_ancombc$run_ancombc(tables[["raw_clade"]],
#                        treatment_key,
#                        dataset_name,
#                        "G",
#                        da_ancombc_dir)

# Indicator Taxa Analysis--------------------------------------------------
# indicsp$run_indicator_analysis(tables[["raw_clade"]],
#                                treatment_key,
#                                dataset_name,
#                                ind_sp_dir)
# indicsp$run_indicator_analysis(tables[["raw_clade"]],
#                                treatment_key,
#                                dataset_name,
#                                ind_sp_dir,
#                                "S")
# indicsp$run_indicator_analysis(tables[["raw_clade"]],
#                                treatment_key,
#                                dataset_name,
#                                ind_sp_dir,
#                                "G")

# Relative Abundance of Individual Crops------------------------------------------------------
 #'*I think tables is more or less the same format as ct already is, so perhaps can just sub in ct instead of tables?*
 #''*don't need "additional taxa", so taking that out (was originally the 4th argument in make_interest_abundance*
plot1 = exploratory$make_interest_abundance(ct,
                                    treatment_key,
                                    dataset_name,
                                    rel_abund_dir) #additional_taxa,

# Relative Abundance of Individual Provinces------------------------------------------------------
exploratory$prov_make_interest_abundance(merged_crops,
                                            treatment_key,
                                            dataset_name,
                                            rel_abund_dir)

#'*requested in meeting*
# Relative Abundance Across All Crops------------------------------------------------------
exploratory$across_crop_abundance(merged_crops,
                                         treatment_key,
                                         dataset_name,
                                         rel_abund_dir)
# Relative Abundance Across All Crops, for one graph for each year (WIP)------------------------------------------------------
exploratory$across_crop_abundance_20(merged_crops,
                                  treatment_key,
                                  dataset_name,
                                  rel_abund_dir)
#' #'*need to group ARGs with small percentages together.  Can adapt the add other bac function? Needs to be after calc_prop*
#' #'*calc_prop produces the .csv file that is saved to the results folder (ie, prop_data* 
#'    # add_other_bac <- function(d) {
#'       column_index <- !(names(ct) %in% c("taxRank","taxLineage","taxID","depth","name"))
#'       # Get the sum counts for taxa of interest, except Bacteria (domain), then
#'       # subtract their sum from Bacteria clade count. This leaves the Bacteria row
#'       # representing all taxa of non-interest
#'       df <- d[d$name!="Bacteria", column_index]
#'       d[d$name=="Bacteria", column_index] <-
#'         d[d$name=="Bacteria", column_index] - as.list(colSums(df, na.rm=T))
#'       # Change the grouping name
#'       d[d$name=="Bacteria",]$name <- "Other ARGs"
#'       
#'       return(d)
#'    # } 
#'     
#'   plot_interest_abundance <- function(d) {
#'     abundance_plot <- ggplot(d, 
#'                              aes(x = crop,
#'                                  y = value,
#'                                  fill = class)) +
#'       geom_bar(stat = "identity", colour = "black") +
#'       facet_grid(~replicate) +
#'       labs(title = "Relative Abundance",
#'            x = "Crop",
#'            y = "Percentage (%)",
#'            fill = "ARG Class") +
#'       theme(axis.text.x = element_text(angle = 45,
#'                                        vjust = 1,
#'                                        hjust = 1))
#'   }

# %>%
#   apply(MARGIN = 1,
#         FUN = function(x) x / sum(x) * 100) %>%
#   t() %>%
#   as.data.frame() %>%
#   mutate(sample_col) %>%
#   relocate(sample)

# Alpha Diversity ---------------------------------------------------------
# exploratory$make_all_alpha_plots(ct,
#                                  treatment_key,
#                                  dataset_name,
#                                  alpha_div_dir)


# # Beta Diversity ----------------------------------------------------------
# exploratory$make_nmds_plots(ct,
#                             treatment_key,
#                             dataset_name,
#                             nmds_dir)
# # _________________________________________________________________________

