# Author(s): Lance Lansing
# Fall 2021
# Normalization via Cumulative Sum Scaling of read count data



# ---------------------------------- Imports -----------------------------------
#import('ggplot2')
import('stringr')
#import('metagenomeSeq')
import('dplyr')
#import('stats', 'setNames', 'aggregate')
# ______________________________________________________________________________

# ------------------------------ Scope variables -------------------------------
firm4 <- paste('Lactobacillus', sep=' ', c('mellifer',
                                           'mellis'))
firm5 <- paste('Lactobacillus', sep=' ', c('kimbladii',
                                           'kullabergensis',
                                           'melliventris',
                                           'helsingborgensis',
                                           'apis')) 
# ______________________________________________________________________________

# --------------------------------- Functions ----------------------------------
export("format_count_table")
format_count_table <- function(tb){
  # Performs some preliminary formatting on count table
  
  # Remove the "Max" and "taxID" columns
  tb <- tb[, !(colnames(tb) %in% c("Max", "Max.1", "taxID"))]
  
  # Taxa name formatting
  tb$name <- gsub(".*<wbr>(.*)", "\\1", tb$name)
  
  # Lineage formatting
  tb$lineage <- gsub("&nbsp;", " ", tb$lineage)
  tb <-  tb %>%
    mutate(lineage = case_when(
      name == 'cellular organisms' ~ 'cellular organisms',
      TRUE ~ str_c(lineage, name, sep = ">")))
  
  # Remove clade read columns 
  tb %>% select(-contains("cladeReads"))
}


export('filter_table')
filter_table <- function(tb){
  # Get only bacterial taxa
  tb <- filter(tb, str_detect(lineage, 'Bacteria'))
}
