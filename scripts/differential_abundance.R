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