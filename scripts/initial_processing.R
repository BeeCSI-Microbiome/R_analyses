# Author(s): Lance Lansing
# Fall 2021
# Normalization via Cumulative Sum Scaling of read count data


# ---------------------------------- Imports -----------------------------------
import('stringr')
import('dplyr')
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
  
  tb <-  tb %>%
    mutate(taxLineage = str_replace(taxLineage, "cellular organisms(>)*", ""))
  
  # Remove clade read columns 
  tb %>% select(-contains("cladeReads"))
}


export('filter_table')
filter_table <- function(tb){
  # Get only bacterial taxa
  tb <- filter(tb, str_detect(taxLineage, 'Bacteria'))
}


export('group_taxa_of_interest')
group_taxa_of_interest <- function(tb){
  tb <- aggregate_unclassified(tb)
  tb <- group_lactobacillus(tb)
}


aggregate_unclassified <- function(tb){
  # Drops 'unclassified <Genus>' taxa, aggregates their counts to the appropriate
  # genus and modifies subtaxa lineage strings accordingly
  
  # Get list of names
  unclass_names <- filter(tb, str_detect(name, 'unclassified'))$name
  genera_names <- str_extract(unclass_names, '(?<=unclassified ).*')
  
  count_cols <- grepl('taxonReads', names(tb))
  
  for (nm in genera_names){
    unclass_nm <- str_c('unclassified ', nm)
    # Get counts of 'unclassified Genus' and 'Genus'
    df <- rbind(tb[tb$name==nm, count_cols],
                tb[tb$name==unclass_nm, count_cols])
    
    # Set counts for Genus to the sum of 'Genus' and 'unclassified Genus'
    tb[tb$name==nm, count_cols] <- as.list(colSums(df, na.rm=T))
    
    # Update lineages of previous 'unclassified Genus' subtaxa
    tb <-  tb %>%
      mutate(taxLineage = case_when(
        str_detect(taxLineage, unclass_nm) ~ gsub(str_c(unclass_nm, '>'),
                                               '', taxLineage),
        
        TRUE ~ taxLineage))
    
    # Drop 'unclassified Genus' row
    tb <- filter(tb, !name==unclass_nm)
  }
  tb
}


group_lactobacillus <- function(tb){
  # Create the Firm-4, Firm-5 and Other Lactobacillus groupings and update the 
  # lineage strings of member species
  firm4_name <- 'Lactobacillus Firm-4'
  firm5_name <- 'Lactobacillus Firm-5'
  other_name <- 'Other Lactobacillus'
  lacto_str <- '>Lactobacillus>'
  lacto_lineage <- tb[tb$name=='Lactobacillus',]$taxLineage
  
  # Return table if this function has already been ran
  if(firm4_name %in% tb$name & firm5_name %in% tb$name & other_name %in% tb$name) {
    return(tb)
  }
  
  # Create rows for Firm-4, Firm-5, and Other
  if(!firm4_name %in% tb$name) {
    tb <- tb %>% add_row(name=firm4_name,
                         taxRank='-',
                         taxLineage=paste(lacto_lineage, firm4_name, sep='>'))
  }
  if(!firm5_name %in% tb$name) {
    tb <- tb %>% add_row(name=firm5_name,
                         taxRank='-',
                         taxLineage=paste(lacto_lineage, firm5_name, sep='>'))
  }
  if(!other_name %in% tb$name) {
    tb <- tb %>% add_row(name=other_name,
                         taxRank='-',
                         taxLineage=paste(lacto_lineage, other_name, sep='>'))
  }
  
  # Update lineage strings of members
  tb <-  tb %>% mutate(taxLineage = case_when(
    name %in% firm4_list ~ 
      gsub(lacto_str, str_c(lacto_str, firm4_name, ">"), taxLineage),
    name %in% firm5_list ~ 
      gsub(lacto_str, str_c(lacto_str, firm5_name, ">"), taxLineage),
    # Other Lactobacillus
    str_detect(taxLineage, lacto_str) & taxRank=='S' ~ 
      gsub(lacto_str, str_c(lacto_str, other_name, ">"), taxLineage),
    TRUE ~ taxLineage))
  
  tb <- reassign_wkb(tb)
}


reassign_wkb <- function(tb){
  # Place 'L. wkB8' within L. helsingborgensis and
  #       'L. wkB10' within L. kullabergensis
  wkb8_str <- 'Lactobacillus sp. wkB8'
  wkb10_str <- 'Lactobacillus sp. wkB10'
  helsing_lineage <- tb[tb$name=='Lactobacillus helsingborgensis',]$taxLineage
  kulla_lineage <- tb[tb$name=='Lactobacillus kullabergensis',]$taxLineage
  
  # Update taxRanks
  tb[tb$name==wkb8_str,]$taxRank <- '-'
  tb[tb$name==wkb10_str,]$taxRank <- '-'
  
  # Correct the lineage strings
  tb[tb$name==wkb8_str,]$taxLineage <- paste(helsing_lineage, wkb8_str, sep='>')
  tb[tb$name==wkb10_str,]$taxLineage <- paste(kulla_lineage, wkb10_str, sep='>')
  tb
}