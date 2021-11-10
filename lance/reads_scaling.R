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
# Reaggregate clade counts with normalized values

# Count lineage_depths. Clade counts will be counted from leaves up
cr_norm <- cr_norm %>%
  mutate(lineage_depth=str_count(lineage, '>'))

# Set remaining NA values to 0
cr_norm[is.na(cr_norm)] <- 0

# Pivot samples into 1 column
cr_norm <- cr_norm %>% 
  tidyr::pivot_longer(
    cols = contains('taxonReads'),
    names_to = "sample",
    values_to = "scaled_reads"
  )

# Clean sample names
cr_norm$sample <- gsub('(.*)\\.taxonReads.y', '\\1', cr_norm$sample)

# Sort rows by lineage depth
cr_norm <- arrange(cr_norm, desc(lineage_depth))

# Sum reads from rows with the following conditions:
#   - samples match
#   - lineage depth is equal (to get taxon count of that taxa), or 1 greater, 
#     to sum up subtaxa
#   - lineage is contained within the subtaxa lineage   
# TODO: This implementation is likely exponential in time. There is probably a better way
for(i in 1:nrow(cr_norm)){
  row <- cr_norm[i,]
  row$scaled_reads <- sum(
    filter(cr_norm,
           sample==row$sample,
           lineage_depth==row$lineage_depth+1 | lineage_depth==row$lineage_depth,
           str_detect(lineage, row$lineage))$scaled_reads)
  cr_norm[i,] <- row
}

# Pivot samples back into their own columns
cr_norm <- cr_norm %>% 
  tidyr::pivot_wider(
    names_from = "sample",
    values_from = c("scaled_reads")
  )

# ______________________________________________________________________________
