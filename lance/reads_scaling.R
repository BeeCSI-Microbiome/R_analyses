# Author(s): Lance Lansing
# Oct 2021
# Normalization via Cumulative Sum Scaling of read count data
# Input: table from Pavian with clades AND taxon counts, not collapsed


# ------------------------------- Package Setup  -------------------------------
packs <- c('ggplot2', 'tidyverse', 'metagenomeSeq', 'plyr')
lapply(packs, library, character.only=TRUE)
# ______________________________________________________________________________


# ---------------------------------- Globals -----------------------------------
# Core list (do not change)
filter_list <- c('root', 'unclassified')

# append any additional taxa (e.g. Apis mellifera)
filter_list <- append(filter_list, c('Apis mellifera'))

css_percentile = 0.5
# ______________________________________________________________________________


# ------------------------------ Table Formatting ------------------------------
# Read in master table
cr <- read_tsv("data/all_clade_and_taxon_reads.tsv")

# Remove the "Max" columns
cr <- cr[, !(colnames(cr) %in% c("Max", "Max.1"))]

# name formatting
cr$name <- gsub(".*<wbr>(.*)", "\\1", cr$name)

# lineage formatting
cr$lineage <- gsub("&nbsp;", " ", cr$lineage)
cr <-  cr %>%
  mutate(lineage = case_when(name=='cellular organisms' ~ 'cellular organisms',
                             TRUE ~ str_c(lineage, name, sep=">")))
# ______________________________________________________________________________


# ------------------------------ Table Filtering -------------------------------
# filter unwanted taxa out
cr <- filter(cr, !name %in% filter_list)
# additional filter needed for uncollapsed clades
# TODO: is there a better way to do this?
if('Apis mellifera' %in% filter_list){
  cr <- filter(cr, !str_detect(cr$lineage, 'Metazoa'))
  }
# ______________________________________________________________________________


# -------------------------------- Split table ---------------------------------
# Create table with only clade data - this will be the unnormalized data
cr_unnorm <- cr %>% select(-contains("taxonReads"))

# Create table with only taxon data but including NA rows for uncollapsed taxa
cr_norm <- cr %>% select(-contains("cladeReads"))
# ______________________________________________________________________________


# ------------------------------- Normalization --------------------------------
# get table with only counts
counts_only <- as.data.frame(select(cr_norm, ends_with('taxonReads')))

# Get table of raw counts with only the taxa that are not ALL NA
# (ie the taxa with clade counts but no taxon count)
# This table will be used in sum visualization later
taxon_raw <- cr_norm[as.logical((rowSums(is.na(counts_only))-ncol(counts_only))),]

# set row names
row.names(counts_only) <- cr_norm$lineage
# remove rows with all NA
counts_only <- counts_only[as.logical((rowSums(is.na(counts_only))-ncol(counts_only))),]
# Set remaining NA values to 0
counts_only[is.na(counts_only)] <- 0

# Perform CSS normalization
css_experiment <- cumNorm(newMRexperiment(counts_only), p=css_percentile)

# Get scaled counts into table
normd <- data.frame(MRcounts(css_experiment, norm=TRUE))
normd <- tibble::rownames_to_column(normd, var='lineage')

# Merge normd back with cr_norm, overwriting count values
cr_norm <- merge(cr_norm, normd, by='lineage', all=T)
# Remove the unnormalized columns
cr_norm <- select(cr_norm, !ends_with('taxonReads.x'))
# ______________________________________________________________________________

# ------------------------------- Visualization --------------------------------
# Call another script which writes a visualization to 'results/scaling_visualization.png'
source('visualize_scaling.r')
# ______________________________________________________________________________


# -------------------------- Calculate Clade Counts --------------------------
# TODO: reaggregate clade counts with normalized values
cr_norm <- cr_norm %>%
  mutate(lineage_depth=str_count(lineage, '>'))
# ______________________________________________________________________________
