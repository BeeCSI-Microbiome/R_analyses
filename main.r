# Author(s): Lance Lansing, Jonathan Ho
# Nov 2021
# Main script for analysis of taxonomic classification results
#
# Input: table from Pavian with clades AND taxon counts, not collapsed


# ------------------------------- Package Setup --------------------------------
packages <- c('tidyverse', 'vegan', 'modules')
lapply(packages, library, character.only = TRUE)
# ______________________________________________________________________________

# ------------------------- Load Aux Files as Modules --------------------------
ip <- use('scripts/initial_processing.R')
scaling <- use('scripts/reads_scaling.R')
exploratory <- use('scripts/exploratory_functions.R')
# ______________________________________________________________________________


# ---------------------------------- Globals -----------------------------------
# Input file paths:
# Counts table (clade and taxon counts, uncollapsed)
counts_path <- "../../2020_ctx_kraken2/ctx_all_clade_taxa_reads_uncollapsed.tsv"
# Percentage table (clade, uncollapsed)
percents_path <- "../../2020_ctx_kraken2/ctx_kraken_all_percent_uncollapsed.tsv" 
# Treatment and replicate names
treat_names <- c("Control", "CLO", "THI")
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")
# Percentile value used by CSS (default=0.5)
css_percentile = 0.5
# ______________________________________________________________________________


# -------------------------------- Read Input ----------------------------------
# Read in count table
ct <- read_tsv(counts_path)
pt <- read_tsv(percents_path)
# ______________________________________________________________________________


# -------------------- Formatting, Filtering, and Scaling ----------------------
# Format and perform filtering on the table
ct <- ip$format_count_table(ct) %>%
  ip$filter_table() %>%
  ip$group_taxa_of_interest()

tables <- scaling$scaling_procedure(ct, css_percentile)

mrexp <- tables$css_MRexp
rt <- tables$raw_taxon
rc <- tables$raw_clade
st <- tables$scaled_taxon
sc <- tables$scaled_clade
# ______________________________________________________________________________


# Relative Abundance ------------------------------------------------------
exploratory$make_interest_abundance(tables[["raw_clade"]],
                                    treat_names,
                                    rep_names)

# Old relative abundance, probably don't need anymore
# exploratory$make_genera_abundance(tables[["raw_clade"]],
#                                   treat_names,
#                                   rep_names)

# Do not use, only for CTX experiment
# exploratory$make_separate_ctx_bars(tables[["raw_clade"]],
#                                    treat_names,
#                                    rep_names)

# Alpha Diversity ---------------------------------------------------------
exploratory$make_all_alpha_plots(tables[["raw_clade"]],
                                 treat_names,
                                 rep_names)


# Beta Diversity ----------------------------------------------------------
exploratory$make_nmds_plots(tables[["raw_clade"]],
                            treat_names,
                            rep_names)

