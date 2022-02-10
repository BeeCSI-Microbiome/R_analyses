<<<<<<< Updated upstream
# Author(s): Lance Lansing
# Jan 2022
# Differential abundance analysis


# ---------------------------------- Imports -----------------------------------
import('data.table')
import('metagenomeSeq')
# ______________________________________________________________________________


# --------------------------------- Functions ----------------------------------
# Adapted from function found in AMR++ shiny app code
export("CSS_normalize_and_extract")
CSS_normalize_and_extract <- function(kraken_MRexp) {
  kraken_norm <- MRcounts(kraken_MRexp, norm=T)
  kraken_raw <- MRcounts(kraken_MRexp, norm=F)
  return(list(kraken_MRexp, kraken_norm, kraken_raw))
}

export('aggregate_and_filter')
aggregate_and_filter <- function(dt_list,
                                 metadata,
                                 filter=TRUE) {
  kraken <- dt_list[[1]]
  kraken_norm <- dt_list[[2]]
  kraken_raw <- dt_list[[3]]
  # <- dt_list[[4]]
  #amr_norm <- dt_list[[5]]
  #amr_raw <- dt_list[[6]]
  #annotations <- dt_list[[7]]
  
  sample_column_id <- colnames(metadata)[1]
  
  # Outputs
  #AMR_analytic_data <- NA
  #AMR_raw_analytic_data <- NA
  #amr_melted_analytic <- NA
  #amr_melted_raw_analytic <- NA
  kraken_analytic_data <- NA
  kraken_raw_analytic_data <- NA
  kraken_melted_analytic <- NA
  kraken_melted_raw_analytic <- NA
  #AMR_analytic_names <- NA
  kraken_analytic_names <- NA
  
  if(!is.na(kraken_norm)) {
    kraken_taxonomy <- data.table(id=rownames(kraken))
    num_sep = max(sapply(kraken_taxonomy$id, function(x) {sum(unlist(strsplit(x, '')) == '>')}))
    
    if(num_sep <= 6) {
      setDT(kraken_taxonomy)[, c('Domain',
                                 'Phylum',
                                 'Class',
                                 'Order',
                                 'Family',
                                 'Genus',
                                 'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
    } else if(num_sep == 7) {
      setDT(kraken_taxonomy)[, c('Domain',
                                 'Kingdom',
                                 'Phylum',
                                 'Class',
                                 'Order',
                                 'Family',
                                 'Genus',
                                 'Species') := tstrsplit(id, '|', type.convert = TRUE, fixed = TRUE)]
    }
    
    kraken_taxonomy[, Kingdom:=NULL]
    
    setkey(kraken_taxonomy, id)
    setkey(kraken_norm, id)
    kraken_norm <- kraken_taxonomy[kraken_norm]  # left outer join
    
    setkey(kraken_raw, id)
    kraken_raw <- kraken_taxonomy[kraken_raw]  # left outer join
    
    
    # Group the kraken data by level for analysis, removing NA entries
    kraken_domain <- kraken_norm[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
    kraken_domain_analytic <- newMRexperiment(counts=kraken_domain[, .SD, .SDcols=!'Domain'])
    rownames(kraken_domain_analytic) <- kraken_domain$Domain
    
    kraken_domain_raw <- kraken_raw[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:8]
    kraken_domain_raw_analytic <- newMRexperiment(counts=kraken_domain_raw[, .SD, .SDcols=!'Domain'])
    rownames(kraken_domain_raw_analytic) <- kraken_domain_raw$Domain
    
    kraken_phylum <- kraken_norm[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
    kraken_phylum_analytic <- newMRexperiment(counts=kraken_phylum[, .SD, .SDcols=!'Phylum'])
    rownames(kraken_phylum_analytic) <- kraken_phylum$Phylum
    
    kraken_phylum_raw <- kraken_raw[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:8]
    kraken_phylum_raw_analytic <- newMRexperiment(counts=kraken_phylum_raw[, .SD, .SDcols=!'Phylum'])
    rownames(kraken_phylum_raw_analytic) <- kraken_phylum_raw$Phylum
    
    kraken_class <- kraken_norm[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
    kraken_class_analytic <- newMRexperiment(counts=kraken_class[, .SD, .SDcols=!'Class'])
    rownames(kraken_class_analytic) <- kraken_class$Class
    
    kraken_class_raw <- kraken_raw[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:8]
    kraken_class_raw_analytic <- newMRexperiment(counts=kraken_class_raw[, .SD, .SDcols=!'Class'])
    rownames(kraken_class_raw_analytic) <- kraken_class_raw$Class
    
    kraken_order <- kraken_norm[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
    kraken_order_analytic <- newMRexperiment(counts=kraken_order[, .SD, .SDcols=!'Order'])
    rownames(kraken_order_analytic) <- kraken_order$Order
    
    kraken_order_raw <- kraken_raw[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:8]
    kraken_order_raw_analytic <- newMRexperiment(counts=kraken_order_raw[, .SD, .SDcols=!'Order'])
    rownames(kraken_order_raw_analytic) <- kraken_order_raw$Order
    
    kraken_family <- kraken_norm[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
    kraken_family_analytic <- newMRexperiment(counts=kraken_family[, .SD, .SDcols=!'Family'])
    rownames(kraken_family_analytic) <- kraken_family$Family
    
    kraken_family_raw <- kraken_raw[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:8]
    kraken_family_raw_analytic <- newMRexperiment(counts=kraken_family_raw[, .SD, .SDcols=!'Family'])
    rownames(kraken_family_raw_analytic) <- kraken_family_raw$Family
    
    kraken_genus <- kraken_norm[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
    kraken_genus_analytic <- newMRexperiment(counts=kraken_genus[, .SD, .SDcols=!'Genus'])
    rownames(kraken_genus_analytic) <- kraken_genus$Genus
    
    kraken_genus_raw <- kraken_raw[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:8]
    kraken_genus_raw_analytic <- newMRexperiment(counts=kraken_genus_raw[, .SD, .SDcols=!'Genus'])
    rownames(kraken_genus_raw_analytic) <- kraken_genus_raw$Genus
    
    kraken_species <- kraken_norm[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
    kraken_species_analytic <- newMRexperiment(counts=kraken_species[, .SD, .SDcols=!'Species'])
    rownames(kraken_species_analytic) <- kraken_species$Species
    
    kraken_species_raw <- kraken_raw[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
    kraken_species_raw_analytic <- newMRexperiment(counts=kraken_species_raw[, .SD, .SDcols=!'Species'])
    rownames(kraken_species_raw_analytic) <- kraken_species_raw$Species
    
    
    # Make long data frame for plotting with ggplot2
    kraken_melted_analytic <- rbind(melt_dt(MRcounts(kraken_domain_analytic), 'Domain'),
                                    melt_dt(MRcounts(kraken_phylum_analytic), 'Phylum'),
                                    melt_dt(MRcounts(kraken_class_analytic), 'Class'),
                                    melt_dt(MRcounts(kraken_order_analytic), 'Order'),
                                    melt_dt(MRcounts(kraken_family_analytic), 'Family'),
                                    melt_dt(MRcounts(kraken_genus_analytic), 'Genus'),
                                    melt_dt(MRcounts(kraken_species_analytic), 'Species'))
    kraken_melted_raw_analytic <- rbind(melt_dt(MRcounts(kraken_domain_raw_analytic), 'Domain'),
                                        melt_dt(MRcounts(kraken_phylum_raw_analytic), 'Phylum'),
                                        melt_dt(MRcounts(kraken_class_raw_analytic), 'Class'),
                                        melt_dt(MRcounts(kraken_order_raw_analytic), 'Order'),
                                        melt_dt(MRcounts(kraken_family_raw_analytic), 'Family'),
                                        melt_dt(MRcounts(kraken_genus_raw_analytic), 'Genus'),
                                        melt_dt(MRcounts(kraken_species_raw_analytic), 'Species'))
    
    # Vector of objects for iteration and their names
    kraken_analytic_data <- c(kraken_domain_analytic,
                              kraken_phylum_analytic,
                              kraken_class_analytic,
                              kraken_order_analytic,
                              kraken_family_analytic,
                              kraken_genus_analytic,
                              kraken_species_analytic)
    kraken_analytic_names <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
    kraken_raw_analytic_data <- c(kraken_domain_raw_analytic,
                                  kraken_phylum_raw_analytic,
                                  kraken_class_raw_analytic,
                                  kraken_order_raw_analytic,
                                  kraken_family_raw_analytic,
                                  kraken_genus_raw_analytic,
                                  kraken_species_raw_analytic)
    
    for( l in 1:length(kraken_analytic_data) ) {
      sample_idx <- match(colnames(MRcounts(kraken_analytic_data[[l]])), metadata[[sample_column_id]])
      pData(kraken_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
      rownames(pData(kraken_analytic_data[[l]])) <- metadata[[sample_column_id]][sample_idx]
      fData(kraken_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_analytic_data[[l]])))
      rownames(fData(kraken_analytic_data[[l]])) <- rownames(MRcounts(kraken_analytic_data[[l]]))
    }
    
    for( l in 1:length(kraken_raw_analytic_data) ) {
      sample_idx <- match(colnames(MRcounts(kraken_raw_analytic_data[[l]])), metadata[[sample_column_id]])
      pData(kraken_raw_analytic_data[[l]]) <- data.frame(
        metadata[sample_idx, .SD, .SDcols=!sample_column_id])
      rownames(pData(kraken_raw_analytic_data[[l]])) <- metadata[sample_idx, .SD, .SDcols=sample_column_id][[sample_column_id]]
      fData(kraken_raw_analytic_data[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_raw_analytic_data[[l]])))
      rownames(fData(kraken_raw_analytic_data[[l]])) <- rownames(MRcounts(kraken_raw_analytic_data[[l]]))
    }
    
    # Ensure that the metadata entries match the factor order of the MRexperiments
    metadata <- data.table(metadata[match(colnames(MRcounts(kraken_domain_analytic)), metadata[[sample_column_id]])])
    setkeyv(metadata, sample_column_id)
    
    # write.csv(make_sparse(kraken_domain, 'Domain', c('Domain')),
    #           'kraken_matrices/sparse_normalized/kraken_Domain_Sparse_Normalized.csv',
    #           row.names=T)
    # write.table(kraken_domain, 'kraken_matrices/normalized/kraken_Domain_Normalized.csv', sep=',', row.names=F, col.names=T)
    # write.table(kraken_domain_raw, 'kraken_matrices/raw/kraken_Domain_Raw.csv', sep=',', row.names=F, col.names=T)
    # 
    # write.csv(make_sparse(kraken_phylum, 'Phylum', c('Phylum')),
    #           'kraken_matrices/sparse_normalized/kraken_Phylum_Sparse_Normalized.csv',
    #           row.names=T)
    # write.table(kraken_phylum, 'kraken_matrices/normalized/kraken_Phylum_Normalized.csv', sep=',', row.names=F, col.names=T)
    # write.table(kraken_phylum_raw, 'kraken_matrices/raw/kraken_Phylum_Raw.csv', sep=',', row.names=F, col.names=T)
    # 
    # write.csv(make_sparse(kraken_class, 'Class', c('Class')),
    #           'kraken_matrices/sparse_normalized/kraken_Class_Sparse_Normalized.csv',
    #           row.names=T)
    # write.table(kraken_class, 'kraken_matrices/normalized/kraken_Class_Normalized.csv', sep=',', row.names=F, col.names=T)
    # write.table(kraken_class_raw, 'kraken_matrices/raw/kraken_Class_Raw.csv', sep=',', row.names=F, col.names=T)
    # 
    # write.csv(make_sparse(kraken_order, 'Order', c('Order')),
    #           'kraken_matrices/sparse_normalized/kraken_Order_Sparse_Normalized.csv',
    #           row.names=T)
    # write.table(kraken_order, 'kraken_matrices/normalized/kraken_Order_Normalized.csv', sep=',', row.names=F, col.names=T)
    # write.table(kraken_order_raw, 'kraken_matrices/raw/kraken_Order_Raw.csv', sep=',', row.names=F, col.names=T)
    # 
    # write.csv(make_sparse(kraken_family, 'Family', c('Family')),
    #           'kraken_matrices/sparse_normalized/kraken_Family_Sparse_Normalized.csv',
    #           row.names=T)
    # write.table(kraken_family, 'kraken_matrices/normalized/kraken_Family_Normalized.csv', sep=',', row.names=F, col.names=T)
    # write.table(kraken_family_raw, 'kraken_matrices/raw/kraken_Family_Raw.csv', sep=',', row.names=F, col.names=T)
    # 
    # write.csv(make_sparse(kraken_genus, 'Genus', c('Genus')),
    #           'kraken_matrices/sparse_normalized/kraken_Genus_Sparse_Normalized.csv',
    #           row.names=T)
    # write.table(kraken_genus, 'kraken_matrices/normalized/kraken_Genus_Normalized.csv', sep=',', row.names=F, col.names=T)
    # write.table(kraken_genus_raw, 'kraken_matrices/raw/kraken_Genus_Raw.csv', sep=',', row.names=F, col.names=T)
    # 
    # write.csv(make_sparse(kraken_species, 'Species', c('Species')),
    #           'kraken_matrices/sparse_normalized/kraken_Species_Sparse_Normalized.csv',
    #           row.names=T)
    # write.table(kraken_species, 'kraken_matrices/normalized/kraken_Species_Normalized.csv', sep=',', row.names=F, col.names=T)
    # write.table(kraken_species_raw, 'kraken_matrices/raw/kraken_Species_Raw.csv', sep=',', row.names=F, col.names=T)
  }
  
  
  
  return(list(AMR_analytic_data,
              AMR_raw_analytic_data,
              amr_melted_analytic,
              amr_melted_raw_analytic,
              AMR_analytic_names,
              kraken_analytic_data,
              kraken_raw_analytic_data,
              kraken_melted_analytic,
              kraken_melted_raw_analytic,
              kraken_analytic_names,
              metadata))
}
=======
# - Differential Abundance --------------------------------------------------

setwd("~/beecsi/R_analyses")

# ------------------------------- Package Setup --------------------------------
packages <- c('tidyverse', 'vegan', 'modules', "data.table", "metagenomeSeq",
              "ggplot2", "here", "PMCMRplus", "broom", "statmod")
lapply(packages, library, character.only = TRUE)
# Adapting one-health continuum code 
# https://github.com/ropolomx/one_health_continuum/blob/master/bcrc_comparative_analysis.R

# Get utility functions
source(here::here('scripts','meg_utility_functions.R'))


# Globals -----------------------------------------------------------------

# File input paths
kraken_analytical <- Sys.glob(here::here("aggregated_data_for_analysis", "cor_2020", "krakenAnalytical_*.csv"))
# Metadata
metadata_filepath = "../data/cor_2020/cor_2020_metadata.csv"
metadata <- read.csv(metadata_filepath, header=T)

# Output paths
# Set the output directory for graphs:
graph_output_dir = 'results/differential_abundance_graphs'
# Set the output directory for statistics:
stats_output_dir = 'results/differential_abundance_stats'

# Create output directories if they don't exist
ifelse(!dir.exists(graph_output_dir), dir.create((graph_output_dir), mode='777'), FALSE)
ifelse(!dir.exists(stats_output_dir), dir.create((stats_output_dir), mode='777'), FALSE)


# Main --------------------------------------------------------------------

kraken_names <- map_chr(
  kraken_analytical,
  ~ str_replace(.x, "^.*_(.*)\\.csv", "\\1")
)

# Read in matrices
temp_kraken_list <- map(
  kraken_analytical,
  ~ read.table(.x, header = T, row.names = 1, sep = ",", quote = "\"")
) %>%
  set_names(nm = kraken_names)

temp_kraken_list <-
  temp_kraken_list %>%
  map(
    ~ mutate(.x, lineage = row.names(.x)) %>% 
      select(.,everything(), lineage) %>%
      gather(key=ID, value = counts, 1:ncol(.x))
  )

temp_kraken_list <- 
  temp_kraken_list %>%
  map(
    ~ group_by(.x, ID, lineage) %>%
      summarise(average_counts = mean(counts))
  )

# Re-widen data-frames
temp_kraken_list <- 
  temp_kraken_list %>%
  map(
    ~ tidyr::spread(.x, key = ID, value = average_counts, fill = 0) 
  )

lineage_row_name <- function(df){
  df <- as.data.frame(df)
  row.names(df) <- df$lineage
  df <- 
    df %>%
    select(-lineage)
  df
}

temp_kraken_list <-
  temp_kraken_list %>%
  map(
    ~ lineage_row_name(.x)
  )

# Create metagenomeSeq MR experiment
kraken_new_mr <-
  temp_kraken_list %>%
  map( ~ newMRexperiment(.x[rowSums(.x) > 0, ]))


# Normalizing unsplit Kraken --------------------------------------

kraken_css <- 
  kraken_new_mr %>%
  map(~ cumNorm(.x,p=0.5))

# Extract the normalized counts into data tables for aggregation

kraken_norm <- 
  kraken_css %>%
  map(~ data.table(MRcounts(.x, norm=T)))

kraken_raw <- 
  kraken_css %>%
  map(~ data.table(MRcounts(.x, norm=F)))

kraken_taxonomy <- 
  temp_kraken_list %>%
  map(~ data.table(id=rownames(.x)))

lineage <- 
  kraken_taxonomy %>%
  map(~ .x$id)

kraken_taxonomy_split <- 
  kraken_taxonomy %>%
  map(~ str_split(string = .x$id, pattern = "\\|"))



# Kraken taxonomy - splitting lineages ----

extract_taxa_of_rank <- function(rank){
  taxa_levels <- list("domain" = "d_",
                      "kingdom" = "k_",
                      "phylum" = "p_",
                      "class" = "c_",
                      "order" = "o_",
                      "family" = "f_",
                      "genus" = "g_",
                      "species" = "s_")
  
  rank_regex <- paste0("^", taxa_levels[rank], ".*")
  
  modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, rank_regex)]) %>%
    modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})
}

species_tax <- modify_depth(kraken_taxonomy_split, .depth = 2, ~ .x[str_detect(.x, "^s_.*")]) %>%
  modify_depth(., .depth = 2, ~ if(length(.x) == 0){.x=NA} else{.x})

domain_tax <- extract_taxa_of_rank("domain")
phylum_tax <- extract_taxa_of_rank("phylum")
class_tax <- extract_taxa_of_rank("class")
order_tax <- extract_taxa_of_rank("order")
family_tax <- extract_taxa_of_rank("family")
genus_tax <- extract_taxa_of_rank("genus")
species_tax <- extract_taxa_of_rank("species")

kraken_tax_dt_clade <- data.table(
  Domain = as.character(domain_tax$cladeReads),
  Phylum = as.character(phylum_tax$cladeReads),
  Class = as.character(class_tax$cladeReads),
  Order = as.character(order_tax$cladeReads),
  Family = as.character(family_tax$cladeReads),
  Genus = as.character(genus_tax$cladeReads),
  Species = as.character(species_tax$cladeReads)
)

kraken_tax_dt_taxon <- data.table(
  Domain = as.character(domain_tax$taxonReads),
  Phylum = as.character(phylum_tax$taxonReads),
  Class = as.character(class_tax$taxonReads),
  Order = as.character(order_tax$taxonReads),
  Family = as.character(family_tax$taxonReads),
  Genus =  as.character(genus_tax$taxonReads),
  Species = as.character(species_tax$taxonReads)
)

kraken_tax_dt <- list(
  "cladeReads" = kraken_tax_dt_clade,
  "taxonReads" = kraken_tax_dt_taxon
)

kraken_tax_dt <- map2(
  kraken_tax_dt,
  lineage,
  ~ .x[, id := .y]
)

kraken_tax_dt <- 
  kraken_tax_dt %>%
  map(
    ~ .x[, lowest := str_split(id, pattern = "\\|") %>% map_chr(~ tail(.x, n=1))]
  )

# Use tax patterns named vector to replace lowest taxon name
# for taxonomy level

tax_levels <- c(
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

tax_regex <- c(
  "^d_.*",
  "^p_.*", 
  "^c_.*", 
  "^o_.*", 
  "^f_.*", 
  "^g_.*", 
  "^s_.*" 
)

# Zip tax levels and tax regex into named vector
tax_patterns <- tax_levels
names(tax_patterns) <- tax_regex

kraken_tax_dt <-
  kraken_tax_dt %>%
  map(
    ~ .x[,lowest_level := str_replace_all(lowest, tax_patterns)]
  )

kraken_tax_dt <- map(kraken_tax_dt, ~ setkey(.x, id))

kraken_norm <- map2(
  kraken_norm,
  kraken_css,
  ~ .x[, id :=(rownames(.y)), ]
)

kraken_norm <- map(
  kraken_norm,
  ~ setkey(.x, id)
)

kraken_norm <- map2(
  kraken_norm,
  kraken_tax_dt,
  ~ .y[.x] # left outer join
)

kraken_norm <-
  kraken_norm %>%
  map( ~ as.data.table(.x))

kraken_raw <- map2(
  kraken_raw,
  kraken_css,
  ~ .x[, id :=(rownames(.y)), ]
)

kraken_raw <- map(
  kraken_raw,
  ~ setkey(.x, id)
)

kraken_raw <- map2(
  kraken_raw,
  kraken_tax_dt,
  ~ .y[.x] # left outer join
) 

kraken_raw <-
  kraken_raw %>%
  map( ~ as.data.table(.x))

# Kraken taxon analytic matrices ------------------------------------------

# Group the kraken taxonReads data by level for analysis, removing NA entries

kraken_taxon_norm_summarised <- 
  tax_levels %>%
  map(
    ~ kraken_norm$taxonReads[!is.na(eval(as.name(.x))) & eval(as.name(.x)) != 'NA' , lapply(.SD, sum), by=.x, .SDcols=!1:10]
  ) %>% 
  set_names(nm=tax_levels)

kraken_taxon_raw_summarised <- 
  tax_levels %>%
  map(
    ~ kraken_raw$taxonReads[!is.na(eval(as.name(.x))) & eval(as.name(.x)) != 'NA', lapply(.SD, sum), by=.x, .SDcols=!1:10]
  ) %>% 
  set_names(nm=tax_levels)

make_analytic <- function(x,y){
  analytic <- newMRexperiment(counts=x[, .SD, .SDcols=!y])
  rownames(analytic) <- x[,eval(as.name(y))]
  analytic
}

kraken_taxon_norm_analytic <- map2(
  kraken_taxon_norm_summarised,
  tax_levels,
  ~ make_analytic(.x, .y)
)

kraken_taxon_raw_analytic <- map2(
  kraken_taxon_raw_summarised,
  tax_levels,
  ~ make_analytic(.x, .y)
)


# Examples for further reference  
# kraken_species_raw <- kraken_raw[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:8]
# kraken_species_raw_analytic <- newMRexperiment(counts=kraken_species_raw[, .SD, .SDcols=!'Species'])
# rownames(kraken_species_raw_analytic) <- kraken_species_raw$Species

# Kraken clade analytic matrices ------------------------------------------

reorder_tax_ranks <- function(level_id){
  level_id <- factor(level_id,
                     levels = c(
                       "Domain",
                       "Phylum",
                       "Class",
                       "Order",
                       "Family",
                       "Genus",
                       "Species"
                     )
  )
  level_id
}



kraken_norm$cladeReads$lowest_level <- reorder_tax_ranks(kraken_norm$cladeReads$lowest_level)
kraken_raw$cladeReads$lowest_level <- reorder_tax_ranks(kraken_raw$cladeReads$lowest_level)
kraken_norm$taxonReads$lowest_level <- reorder_tax_ranks(kraken_norm$taxonReads$lowest_level)
kraken_raw$taxonReads$lowest_level <- reorder_tax_ranks(kraken_raw$taxonReads$lowest_level)

kraken_clade_norm_list <-
  kraken_norm$cladeReads %>%
  split(.$lowest_level) %>%
  map(
    ~ .x %>%
      select(-c(Domain:id,lowest_level))
  )

kraken_clade_raw_list <-
  kraken_raw$cladeReads %>%
  split(.$lowest_level) %>%
  map(
    ~ .x %>%
      select(-c(Domain:id,lowest_level))
  )

make_analytic_clade <- function(x){
  analytic <- newMRexperiment(counts=x[, .SD, .SDcols=!'lowest'])
  rownames(analytic) <- x$lowest
  analytic
}

kraken_clade_norm_analytic <- map(
  kraken_clade_norm_list,
  ~ make_analytic_clade(.x)
)

kraken_clade_raw_analytic <- map(
  kraken_clade_raw_list,
  ~ make_analytic_clade(.x)
)

# Make long data frame for plotting with ggplot2

kraken_taxon_norm_melted <- imap_dfr(
  kraken_taxon_norm_analytic,
  ~ melt_dt(MRcounts(.x), .y)
) # getting warning: binding character and factor vector, coercing into character vector

kraken_taxon_raw_melted <- imap_dfr(
  kraken_taxon_raw_analytic,
  ~ melt_dt(MRcounts(.x), .y)
) # getting warning: binding character and factor vector, coercing into character vector

kraken_clade_norm_melted <- imap_dfr(
  kraken_clade_norm_analytic,
  ~ melt_dt(MRcounts(.x), .y)
)

kraken_clade_raw_melted <- imap_dfr(
  kraken_clade_raw_analytic,
  ~ melt_dt(MRcounts(.x), .y)
) # getting warning: binding character and factor vector, coercing into character vector


# kraken_taxon_norm_melted$Level_ID <- reorder_tax_ranks(kraken_taxon_norm_melted$Level_ID)
# kraken_taxon_raw_melted$Level_ID <- reorder_tax_ranks(kraken_taxon_raw_melted$Level_ID)
# kraken_clade_norm_melted$Level_ID <- reorder_tax_ranks(kraken_clade_norm_melted$Level_ID)
# kraken_clade_raw_melted$Level_ID <- reorder_tax_ranks(kraken_clade_raw_melted$Level_ID)

# Match metadata ----------------------------------------------------------
metadata <- data.table(metadata[match(colnames(MRcounts(kraken_clade_norm_analytic[[1]])), metadata[, "ID"]), ])


# Pair metadata with kraken data
for( l in 1:length(kraken_taxon_norm_analytic) ) {
  sample_idx <- match(colnames(MRcounts(kraken_taxon_norm_analytic[[l]])), metadata[["ID"]])
  pData(kraken_taxon_norm_analytic[[l]]) <- data.frame(
    metadata[sample_idx, .SD, .SDcols=!"ID"])
  rownames(pData(kraken_taxon_norm_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols="ID"][["ID"]]
  fData(kraken_taxon_norm_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_taxon_norm_analytic[[l]])))
  rownames(fData(kraken_taxon_norm_analytic[[l]])) <- rownames(MRcounts(kraken_taxon_norm_analytic[[l]]))
  pData(kraken_taxon_norm_analytic[[l]]@expSummary$expSummary)$normFactors <- calcNormFactors(kraken_new_mr$cladeReads, p=0.5)
}

for( l in 1:length(kraken_taxon_raw_analytic) ) {
  sample_idx <- match(colnames(MRcounts(kraken_taxon_raw_analytic[[l]])), metadata[["ID"]])
  pData(kraken_taxon_raw_analytic[[l]]) <- data.frame(
    metadata[sample_idx, .SD, .SDcols=!"ID"])
  rownames(pData(kraken_taxon_raw_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols="ID"][["ID"]]
  fData(kraken_taxon_raw_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_taxon_raw_analytic[[l]])))
  rownames(fData(kraken_taxon_raw_analytic[[l]])) <- rownames(MRcounts(kraken_taxon_raw_analytic[[l]]))
}

for( l in 1:length(kraken_clade_norm_analytic) ) {
  sample_idx <- match(colnames(MRcounts(kraken_clade_norm_analytic[[l]])), metadata[["ID"]])
  pData(kraken_clade_norm_analytic[[l]]) <- data.frame(
    metadata[sample_idx, .SD, .SDcols=!"ID"])
  rownames(pData(kraken_clade_norm_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols="ID"][["ID"]]
  fData(kraken_clade_norm_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_clade_norm_analytic[[l]])))
  rownames(fData(kraken_clade_norm_analytic[[l]])) <- rownames(MRcounts(kraken_clade_norm_analytic[[l]]))
  pData(kraken_clade_norm_analytic[[l]]@expSummary$expSummary)$normFactors <- calcNormFactors(kraken_new_mr$taxonReads, p=0.5)
}

for( l in 1:length(kraken_clade_raw_analytic) ) {
  sample_idx <- match(colnames(MRcounts(kraken_clade_raw_analytic[[l]])), metadata[["ID"]])
  pData(kraken_clade_raw_analytic[[l]]) <- data.frame(
    metadata[sample_idx, .SD, .SDcols=!"ID"])
  rownames(pData(kraken_clade_raw_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols="ID"][["ID"]]
  fData(kraken_clade_raw_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_clade_raw_analytic[[l]])))
  rownames(fData(kraken_clade_raw_analytic[[l]])) <- rownames(MRcounts(kraken_clade_raw_analytic[[l]]))
}

kraken_taxon_names <- names(kraken_taxon_raw_analytic)
kraken_clade_names <- names(kraken_clade_raw_analytic)


# Statistical Analysis Initialization -------------------------------------

statistical_analyses = list(
  # ACTIVITY 2 - CORN (likely applicable to most activity 2 datasets)
  # Description: , control for replicate using random effect
  list(
    name = 'Exposure',
    subsets = list(),
    model_matrix = '~ 0 + Exposed',
    contrasts = list(
      'ExposedTRUE - ExposedFALSE'
    ),
    random_effect = 'Replicate'
  )
)

# statistical_analyses = list(
#   # ACTIVITY 1 - CTX
#   # Description: , control for replicate using random effect
#   list(
#     name = 'Treatment',
#     subsets = list(),
#     model_matrix = '~ 0 + Treatment',
#     contrasts = list(
#       'Treatmentcontrol - Treatmentclothianidin',
#       'Treatmentcontrol - Treatmentthiamethoxam',
#       'Treatmentclothianidin - Treatmentthiamethoxam'
#     ),
#     random_effect = 'Replicate'
#   )
# )

source(here('scripts','meg_utility_functions.R'))

# Apply differential abundance analysis
for (a in 1:length(statistical_analyses)){
  meg_fitZig(data_list=kraken_taxon_norm_analytic[2:7],
             data_names=kraken_taxon_names[2:7],
             metadata=metadata,
             zero_mod=model.matrix(~1 + log(libSize(kraken_css$taxonReads))),
             data_mod=statistical_analyses[[a]]$model_matrix,
             filter_min_threshold=0.15,
             contrast_list=statistical_analyses[[a]]$contrasts,
             random_effect_var=statistical_analyses[[a]]$random_effect,
             outdir=stats_output_dir,
             analysis_name=statistical_analyses[[a]]$name,
             analysis_subset=statistical_analyses[[a]]$subsets,
             data_type='Microbiome_taxonReads',
             pval=0.1,
             top_hits=1000)
}

for (a in 1:length(statistical_analyses)){
  meg_fitZig(data_list=kraken_clade_norm_analytic[2:7],
             data_names=kraken_clade_names[2:7],
             metadata=metadata,
             zero_mod=model.matrix(~1 + log(libSize(kraken_css$cladeReads))),
             data_mod=statistical_analyses[[a]]$model_matrix,
             filter_min_threshold=0.15,
             contrast_list=statistical_analyses[[a]]$contrasts,
             random_effect_var=statistical_analyses[[a]]$random_effect,
             outdir=stats_output_dir,
             analysis_name=statistical_analyses[[a]]$name,
             analysis_subset=statistical_analyses[[a]]$subsets,
             data_type='Microbiome_cladeReads',
             pval=0.1,
             top_hits=1000)
}
>>>>>>> Stashed changes
