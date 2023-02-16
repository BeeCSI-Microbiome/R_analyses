# This code snippet may be used with the 'raw clade' table produced in main.R
# to produce a list of taxa above a certain abundance threshold.

# The min_abundance is the percent at which a taxon's relative abundance
# must be greater or equal to in at least n samples, where n is the value
# specified as min_samples_above_cutoff
taxa_cutoff_explore <- function(taxa_count_matrix, dataset_name) {
  # Minimum percent abundance for a taxon
  min_abundance <- 1
  # Minimum number of samples that a taxon must exceed the abundance cutoff
  min_samples_above_cutoff <- 2


  pc_sc <- select(taxa_count_matrix, !matches('(taxID|taxLineage|depth)'))
  bac_sc <- pc_sc %>% filter(name == 'Bacteria')

  # Get proportions using Bacteria count as total
  pc_sc[c(-1,-2)] <- sapply(X=colnames(pc_sc[c(-1,-2)]),
                            FUN = function(colname){
                              col <- pc_sc[,colname]
                              #print(bac_sc[[colname]])
                              col / bac_sc[[colname]] * 100
                            })

  taxa_of_interest_regex <- "Paenibacillus larvae|Melissococcus plutonius|Bartonella apis|Snodgrassella alvi|Lactobacillus|Frischella perrara|Bifidobacterium|Gilliamella apicola|Gilliamella apis| sp\\."

  # Get table of taxa above cutoff
  taxa_above_cutoff <- pc_sc %>%
    mutate(num_samples_above_thresh = rowSums(pc_sc[c(-1,-2)] >= min_abundance)) %>%
    filter(num_samples_above_thresh >= min_samples_above_cutoff &
           taxRank == "S") %>%
    mutate(new_taxa_of_interest = !str_detect(name, taxa_of_interest_regex))

  xl_file <- loadWorkbook("additional_taxa_above_1p.xlsx")
  tryCatch({addWorksheet(xl_file, dataset_name)},
           error=function(cond){message(cond)})

  writeData(xl_file, sheet = dataset_name,
            x = select(taxa_above_cutoff, name, taxRank, new_taxa_of_interest))
  saveWorkbook(xl_file, "additional_taxa_above_1p.xlsx", overwrite = TRUE)

  return(taxa_above_cutoff)
}
