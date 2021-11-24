# Author(s): Lance Lansing
# Nov 2021
# Main script for analysis of taxonomic classification results
#
# Input: table from Pavian with clades AND taxon counts, not collapsed


# ------------------------------- Package Setup --------------------------------
#packs <- c('ggplot2', 'tidyverse', 'metagenomeSeq', 'plyr')
packages <- c('tidyverse', 'modules')
lapply(packages, library, character.only = TRUE)
# ______________________________________________________________________________

# ------------------------- Load Aux Files as Modules --------------------------
scaling <- use('scripts/reads_scaling.R')
# ______________________________________________________________________________


# ---------------------------------- Globals -----------------------------------
# Input file paths:
# Counts table (clade and taxon counts, uncollapsed)
counts_path <- "data/all_clade_and_taxon_reads.tsv"
# Percentage table (clade, uncollapsed)
percents_path <- "data/all_clades_percents.tsv" 

# Core list (do not change)
filter_list <- c('root', 'unclassified') %>%
# append any additional taxa (e.g. Apis mellifera)
  append(c('Apis mellifera'))

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
  scaling$filter_table(filter_list)

# Get just the taxon counts by dropping taxa with no counts in any sample
ct_raw_taxon <- scaling$drop_all_NA_rows(ct)
# Get clade counts. We calculate counts rather than use clade counts from Pavian
# in order to account for taxa that we filtered out
ct_raw_clade <- scaling$get_raw_clade_data(ct)

# Perform the scaling
ct_scaled <- scaling$css_scale(ct, css_percentile)

# Get scaled taxon table by dropping taxa rows with no counts in any sample
ct_scaled_taxon <- scaling$drop_all_NA_rows(ct_scaled)

# Calculate scaled clade counts from scaled taxon counts 
ct_scaled_clade <- scaling$calc_clade_counts(ct_scaled)
# ______________________________________________________________________________