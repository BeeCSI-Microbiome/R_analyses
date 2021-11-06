# Author(s): Lance Lansing
# Oct 2021
# Produce a summary of read information

# ----------------------------------- Setup ------------------------------------
packs <- c('ggplot2', 'tidyverse', 'metagenomeSeq', 'plyr')
lapply(packs, library, character.only=TRUE)
# ------------------------------------------------------------------------------

# ---------------------------------- Globals -----------------------------------
# Core list (do not change)
taxa_list <- c('root', 'unclassified','Apis mellifera')

# append any additional taxa (e.g. Lactobacillus)
taxa_list <- append(taxa_list, c('Lactobacillus', 'Gilliamela'))

# Is this needed?
taxa_list_nospace <- gsub(' ', '_', taxa_list)
# ------------------------------------------------------------------------------


# ------------------------------ Table Formatting ------------------------------
# Read in master table
cr <- read_tsv("data/all_clade_and_taxon_reads.tsv")

# Pivot samples to a column, keep only taxa in taxa_list
ta <- subset(cr, name %in% taxa_list) %>% 
  select_if(str_detect(names(.), 'name|cladeReads')) %>% 
  tidyr::pivot_longer(
    cols = ends_with('cladeReads'),
    names_to = "sample",
    values_to = "reads"
  )

# sample name formatting
ta$sample <- gsub('(.*)\\.cladeReads', '\\1', ta$sample)

# Pivot to get taxa as columns
ta <- ta %>% 
  tidyr::pivot_wider(
    names_from = "name",
    values_from = c("reads")
  )

# some column name formatting
ta <- dplyr::rename(ta, classified=root)
ta <- dplyr::rename_with(ta, ~ gsub(' ', '_', .x))
ta <- dplyr::rename_with(ta, ~ str_c(.x, '_reads'), !matches('sample'))

# create total read count column
ta <- dplyr::mutate(ta, total_reads = unclassified_reads+classified_reads, .after=sample)

# create "percent of total reads" columns
ta <- ta %>% dplyr::mutate(across(
  !matches('sample|total_reads'),
  ~ .x/total_reads*100,
  .names = "percent_total_{.col}",
  ))

# create "percent of classified reads" columns
ta <- ta %>% dplyr::mutate(across(
  !matches('sample|total_reads|*classified*|*percent*'),
  ~ .x/classified_reads*100,
  .names = "percent_classified_{.col}",
))

# column name formatting
ta <- dplyr::rename_with(ta, ~ gsub('(percent.*)_reads', '\\1', .x))

write.table(ta, file='results/read_summary.tsv', quote=F, sep='\t', row.names=F)
