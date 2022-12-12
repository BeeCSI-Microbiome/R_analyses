# Normalization via Cumulative Sum Scaling of read count data

# ---------------------------------- Imports -----------------------------------
import('stringr')
import('dplyr')
import('R6')
import('data.tree')
import('treemap')
import('DiagrammeR')
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
#'* it's assigning tb as ct, which is from main.r and is the data matrix (after some reformatting)*
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


export('calculate_clade_counts')
#'*being fed ct*
calculate_clade_counts <- function(tb){
  # Given the formatted and filtered table, return a list of two items:
  # (Raw taxon table, raw clade table)
  
  #'*this is basically removing NA and then removing the taxonReads prefix from sample names, which has already been done for our data*
  tb_raw_taxon <- drop_all_NA_rows(tb)
  names(tb_raw_taxon) <- gsub("taxonReads_", "", names(tb_raw_taxon))
  
  # We calculate counts rather than use clade counts from Pavian in order to 
  # account for taxa that we filtered out
  tb_raw_clade <- calc_clade_counts(tb)
  
  list(raw_taxon=tb_raw_taxon,
       raw_clade=tb_raw_clade)
}


drop_all_NA_rows <- function(tb){
  # drops all rows from the table in which counts are NA for all samples
  counts_only <- as.data.frame(select(tb, contains('taxonReads')))
  tb <- tb[as.logical((rowSums(is.na(counts_only))-ncol(counts_only))), ]
  tb[is.na(tb)] <- 0
  tb
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


calc_clade_counts <- function(tb) {
  
  # Notes can be cleaned up but may reduce clarity. 

  #Pseudo-Code (Rough)
  # 1. Clean data (just reads and tree key)
  # 2. Get string vector for column names
  # 3. Convert data to tree.
  # 4. make dataframe w/ name 
  # e.g. temp <- as.data.frame(tb_tree$Get(attribute = "name"))
  # 5. iterate through aggregate function (#see testing #1-4) # should self clear so probably don't have to reset aggreagate data structure
  # 6. combine Aggregate data frames w/ the variable name. 
  # 7. Return dataframe. 
  
  
  # 1. Cleaning
  counts_only <- as.data.frame(select(tb, contains(c('taxonReads', 'taxLineage', 'depth', "taxRank","taxID"))))
  tb_raw_taxon <- counts_only[as.logical((rowSums(is.na(counts_only))-ncol(counts_only))), ]
  tb_raw_taxon[is.na(tb_raw_taxon)] <- 0
  
  names(tb_raw_taxon) <- gsub("taxonReads_", "", names(tb_raw_taxon))
  
  tb_raw_taxon$pathString <- str_replace_all(tb_raw_taxon$taxLineage, '/', # must be path string
                                             '-')
  
  tb_raw_taxon$pathString <- str_replace_all(tb_raw_taxon$pathString, '>', # must be path string
                                             '/')
  # 2 Column names vector 
  columnNames <- colnames(tb_raw_taxon)
  truncated_colNames <- columnNames[columnNames %in% c("taxRank","taxID", "depth","taxLineage", "pathString") == FALSE]
  
  # 3 convert data to tree
  tb_tree <- as.Node(tb_raw_taxon)
  
  # 4 name data frame <not optimal but hopefully not too costly>
  temp1 <- as.data.frame(tb_tree, row.names = NULL, optional = FALSE, "name") # join on name
  
  temp3 <- temp1 # this is redundant but helps with troubleshooting for now 
  
  # join taxRank
  temp2 <- as.data.frame(tb_tree, row.names = NULL, optional = FALSE, "name", "taxRank") # join on name
  
  temp3 <- merge(temp3, temp2)
  
  #join taxID
  temp2 <- as.data.frame(tb_tree, row.names = NULL, optional = FALSE, "name", "taxID") # join on name
  
  temp3 <- merge(temp3, temp2)
  
  # join depth
  temp2 <- as.data.frame(tb_tree, row.names = NULL, optional = FALSE, "name", "depth") # join on name
  
  temp3 <- merge(temp3, temp2)
  
  # currently broken check merge 
  # 5 / 6 work simultaneously
  for(colName in truncated_colNames) {
    tb_tree$Do(function(node) {
      node$Aggregation <- sum(if (is.na(GetAttribute(node, attribute = colName))) 0
                              else GetAttribute(node, attribute = colName), 
                              if (node$isLeaf) 0 else data.tree::Aggregate(node, attribute = "Aggregation", aggFun = sum)
      )}, traversal = "post-order")
    
    temp2 <- as.data.frame(tb_tree, row.names = NULL, optional = FALSE, "name", "Aggregation") # join on name
    colnames(temp2) <- c("levelName", "name", colName)
    
    temp3 <- merge(temp3, temp2)
  }

  # Joining taxlineage for fully correct re-factoring
  temp2 <- as.data.frame(tb_tree, row.names = NULL, optional = FALSE, "name", "taxLineage") # join on name
  
  temp3 <- merge(temp3, temp2)
  
  # final formatting.
  temp3 <- subset(temp3, select = -levelName)
  
  temp3 <- subset(temp3, name != "root2")
  
  return(temp3)
}
