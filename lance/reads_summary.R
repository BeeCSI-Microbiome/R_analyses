# Author(s): Lance Lansing
# Oct 2021
# Produce a summary of read information


# ----------------------------------- Setup  -----------------------------------
packs <- c('ggplot2', 'tidyverse', 'metagenomeSeq', 'plyr')
lapply(packs, library, character.only=TRUE)
# ------------------------------------------------------------------------------


# ------------------------------ Table Formatting ------------------------------
# Read in master table
cr <- read_tsv("data/all_clade_and_taxon_reads.tsv")

# Pivot samples to a column, keep only 3 taxa seen below
ta <- cr %>% subset(str_detect(name, 'unclassified$|cellular organisms|Apis mellifera')) %>% 
  tidyr::pivot_longer(
    cols = 5:14,
    names_to = "sample",
    values_to = "reads"
  )

# Remove the "Max" columns
ta <- ta[, (colnames(ta) %in% c("sample", "reads", "name"))]
ta$sample <- gsub('(.*)\\.cladeReads', '\\1', ta$sample)

# Pivot to get taxa as columns
ta <- ta %>% 
  tidyr::pivot_wider(
    names_from = "name",
    values_from = c("reads")
  )

# column name formatting
ta <- setNames(ta, c('sample','unclassified_reads','classified_reads','Apis_mellifera_reads'))  

# mutate to create additional columns
# Total and percents
ta <- mutate(ta, total_reads = unclassified_reads+classified_reads) 
ta <- mutate(ta, percent_classified = classified_reads/total_reads*100)
ta <- mutate(ta, percent_unclassified = unclassified_reads/total_reads*100)
ta <- mutate(ta, percent_Apis_mellifera = Apis_mellifera_reads/total_reads*100)
# meta data 
ta <- mutate(ta, exposure=ifelse(endsWith(sample, "e"), "exposed", "unexposed"))
ta <- mutate(ta, replicate=str_extract(sample, '\\d\\d'))

# Reorder columns
ta <- ta[, c(1,9,10,5,2,7,3,6,4,8)]

write.table(ta, file='results/read_summary.tsv', quote=F, sep='\t', row.names=F)
