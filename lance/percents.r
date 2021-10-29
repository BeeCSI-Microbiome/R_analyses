# Author(s): Lance Lansing
# Oct 2021
# This script takes taxonomic % data and produces some visualizations

packs <- c('ggplot2', 'tidyverse')
lapply(packs, library, character.only=TRUE)
# --------- Table Formatting
# Read in master table
acp <- read_tsv("all_clades_percent.tsv")

# Remove the "Max" column
acp <- acp[, !(colnames(acp) %in% c("Max"))]

# pivot longer to get sample as a 
acp <- acp %>% 
  tidyr::pivot_longer(
    cols = 4:13,
    names_to = "sample",
    values_to = "percent"
  )

# Clean up the sample name entries
acp$sample <- gsub('(\\w+\\d\\d(u|e)).*$', "\\1", acp$sample)

# Create the "exposure" column (assumes sample ends with 'u' or 'e')
acp <- mutate(acp, exposure=ifelse(endsWith(sample, "e"), "exposed", "unexposed"))
# Create the "replicate" column
acp <- mutate(acp, replicate=str_extract(sample, '\\d\\d'))

# Get separate taxa level tables
species_percents <- filter(acp, str_detect(taxRank, 'S'))
genus_percents <- filter(acp, str_detect(taxRank, 'G'))
family_percents <- filter(acp, str_detect(taxRank, 'F'))
order_percents <- filter(acp, str_detect(taxRank, 'O'))
class_percents <- filter(acp, str_detect(taxRank, 'C'))
phylum_percents <- filter(acp, str_detect(taxRank, 'P'))
kingdom_percents <- filter(acp, str_detect(taxRank, 'K'))
domain_percents <- filter(acp, str_detect(taxRank, 'D'))
unclassified_percents <- filter(acp, str_detect(taxRank, 'U'))


violin_p <- function(df, taxa_name, x_meta, color_meta=NULL) {
  # This function returns a violin plots separated by the given x_meta column of
  # the given data frame df. It displays % of reads of the given taxa_name.
  # Dots can be colored according to a given column color_meta (optional)
  plot <- subset(df, name==taxa_name) %>%
    ggplot(aes(x={{x_meta}}, y=percent)) +
    geom_violin() +
    #ylim(0, 10) +
    labs(title=str_c('Percent of reads that are ', taxa_name),
         x='Metadata category',
         y=str_c('Percent ', taxa_name)) +
    geom_dotplot(binaxis = 'y',
                 dotsize = 1,
                 stackdir = 'center',
                 stackgroups=T,
                 binpositions='all',
                 aes(fill={{color_meta}}))
  return(plot)
}

unc_vp <- acp %>% 
  violin_p(taxa_name = 'unclassified',
           x_meta = exposure,
           color_meta = replicate)
unc_vp

apis_mellifera_vp <- acp %>% 
  violin_p(taxa_name = 'Apis mellifera',
           x_meta = exposure,
           color_meta = replicate)
apis_mellifera_vp

nosema_vp <- acp %>% 
  violin_p(taxa_name = 'Nosema',
           x_meta = exposure,
           color_meta = replicate)
nosema_vp

