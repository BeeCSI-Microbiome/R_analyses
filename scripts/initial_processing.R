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
firm4_list <- paste('Lactobacillus', sep=' ', c('mellifer',
                                           'mellis'))
firm5_list <- paste('Lactobacillus', sep=' ', c('kimbladii',
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


export('group_lactobacillus')
group_lactobacillus <- function(tb){
  firm4_name <- 'Lactobacillus Firm-4'
  firm5_name <- 'Lactobacillus Firm-5'
  lacto_str <- '>Lactobacillus>'
  lacto_lineage <- tb[tb$name=='Lactobacillus',]$lineage
  # Create rows for Firm-4 and Firm-5 
  tb <- tb %>% add_row(name=firm4_name,
                     taxRank='-',
                     lineage=paste(lacto_lineage, firm4_name, sep='>'))
  tb <- tb %>% add_row(name=firm5_name,
                     taxRank='-',
                     lineage=paste(lacto_lineage, firm5_name, sep='>'))
  
  tb <-  tb %>%
    mutate(lineage = case_when(
      name %in% firm4_list ~ gsub(lacto_str,
                                  str_c(lacto_str, firm4_name, ">"),
                                  lineage),
      name %in% firm5_list ~ gsub(lacto_str,
                                  str_c(lacto_str, firm5_name, ">"),
                                  lineage),
      TRUE ~ lineage))
}
