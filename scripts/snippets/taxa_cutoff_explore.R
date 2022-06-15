# This code snippet may be used with the 'raw clade' table produced in main.R
# to produce a list of taxa above and below a certain abundance threshold.

# The min_abundance is the percent at which a taxon's relative abundance 
# must be greater or equal to in at least n samples, where n is the value 
# specified as min_samples_above_cutoff

# Minimum percent abundance for a taxon
min_abundance <- 1
# Minimum number of samples that a taxon must exceed the abundance cutoff
min_samples_above_cutoff <- 2


pc_sc <- select(tables[["raw_clade"]], matches('.*(Reads|name|Rank).*'))
bac_sc <- pc_sc %>% filter(name == 'Bacteria')

# Get proportions using Bacteria count as total
pc_sc[c(-1,-2)] <- sapply(X=colnames(pc_sc[c(-1,-2)]),
                      FUN = function(colname){
                        col <- pc_sc[,colname]
                        #print(bac_sc[[colname]])
                        mutate(col / bac_sc[[colname]] * 100)
                      })

# Get table of taxa above cutoff
taxa_above_cutoff <- pc_sc %>% 
  mutate(num_samples_above_thresh = rowSums(pc_sc[c(-1,-2)] >= min_abundance)) %>% 
  filter(num_samples_above_thresh >= min_samples_above_cutoff)

# Get table of taxa below cutoff
# taxa_below_cutoff <- pc_sc %>% 
#   mutate(num_samples_above_thresh = rowSums(pc_sc[c(-1,-2)] >= min_abundance)) %>% 
#   filter(num_samples_above_thresh < min_samples_above_cutoff)

write_csv(taxa_above_cutoff[c(1,2)], file=str_glue("{main_outdir}/taxa_above_cutoff_{min_abundance}p.csv"))
# write_csv(taxa_below_cutoff[c(1,2)], file=str_glue("{main_outdir}/taxa_below_cutoff_{min_abundance}p.csv"))
