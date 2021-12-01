# Author(s): Lance Lansing, Jonathan Ho
# Nov 2021
# Main script for analysis of taxonomic classification results
#
# Input: table from Pavian with clades AND taxon counts, not collapsed


# ------------------------------- Package Setup --------------------------------
#packs <- c('ggplot2', 'tidyverse', 'metagenomeSeq', 'plyr')
packages <- c('tidyverse', 'vegan', 'modules')
lapply(packages, library, character.only = TRUE)
# ______________________________________________________________________________

# ------------------------- Load Aux Files as Modules --------------------------
scaling <- use('scripts/reads_scaling.R')
exploratory <- use('scripts/exploratory_functions.R')
# ______________________________________________________________________________


# ---------------------------------- Globals -----------------------------------
# Input file paths:
# Counts table (clade and taxon counts, uncollapsed)
counts_path <- "../2020_ctx_kraken2/ctx_all_clade_taxa_reads_uncollapsed.tsv"
# Percentage table (clade, uncollapsed)
percents_path <- "../2020_ctx_kraken2/ctx_kraken_all_percent_uncollapsed.tsv" 

# Core list (do not change)
#filter_list <- c('root', 'unclassified') %>%
# append any additional taxa (e.g. Apis mellifera)
#  append(c('Apis mellifera'))

css_percentile = 0.5
# ______________________________________________________________________________


# -------------------------------- Read Input ----------------------------------
# Read in count table
ct <- read_tsv(counts_path)
pt <- read_tsv(percents_path)
# ______________________________________________________________________________


# -------------------- Formatting, Filtering, and Scaling ----------------------
# Format and perform filtering on the table
ct <- scaling$format_count_table(ct) %>%
  scaling$filter_table()

tables <- scaling$scaling_procedure(ct, css_percentile)
# ______________________________________________________________________________


# Relative Abundance ------------------------------------------------------
# TODO: these details will need to be provided
treat_names <- c("Control", "CLO", "THI")
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")
plot_title <- "CTX Abundance Using Percent(%) Data"

relative_abundance <- tables[["raw_clade"]] %>%
  exploratory$make_genera_abundance(treat_names,
                                    rep_names,
                                    plot_title)
relative_abundance

# Alpha Diversity ---------------------------------------------------------
# TODO: this will need to be defined
alpha_index <- "Shannon"

alpha_div <- tables[["raw_clade"]] %>%
  exploratory$make_alpha_bars(treat_names,
                              rep_names,
                              alpha_index)
alpha_div
