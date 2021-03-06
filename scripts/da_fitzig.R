# - Differential Abundance --------------------------------------------------

# Adapting one-health continuum code 
# https://github.com/ropolomx/one_health_continuum/blob/master/bcrc_comparative_analysis.R

# Get utility functions

import("glue")
import("data.table")
import("metagenomeSeq")
import("utils")
import("purrr")
import("stringr")
import("tidyr")
import("dplyr")
import("Biobase")
import("stats")
import("tibble")

meg_functions <- use("scripts/meg_utility_functions.R")


# This function is for recalculating clade counts from scaled taxon counts
calc_clade_counts <- function(tb){
  
  # make row names into lineage column
  tb <- cbind(rownames(tb), data.frame(tb, row.names=NULL))
  colnames(tb)[1] <- "lineage"
  
  # Count lineage_depth. Clade counts will be counted from the deepest leaves up
  tb <- tb %>%
    mutate(lineage_depth = str_count(lineage, '\\|'), .after = lineage)
  
  # '|' characters must be changed because they interfere with regex functions
  tb <- tb %>% 
    mutate(lineage = str_replace_all(lineage, "\\|", ">"))
  
  # Set remaining NA values to 0
  tb[is.na(tb)] <- 0
  
  # Pivot samples into 1 column
  tb <- tb %>% 
    tidyr::pivot_longer(
      cols = -c(lineage, lineage_depth),
      names_to = "sample",
      values_to = "scaled_reads"
    )
  
  # Sort rows by lineage depth
  tb <- arrange(tb, desc(lineage_depth))
  
  # Sum reads from rows with the following conditions:
  #   - samples match
  #   - lineage depth is equal (to get taxon count of that taxa), or 1 greater, 
  #     to sum up subtaxa
  #   - lineage is contained within the subtaxa lineage   
  # TODO: This implementation is likely exponential in time. There is probably a better way
  for(i in 1:nrow(tb)){
    row <- tb[i, ]
    row$scaled_reads <- sum(
      dplyr::filter(tb,
             sample == row$sample,
             lineage_depth == row$lineage_depth + 1 | lineage_depth == row$lineage_depth,
             str_detect(lineage, row$lineage))$scaled_reads)
    tb[i, ] <- row
  }
  
  # Pivot samples back into their own columns
  tb <- tb %>% 
    tidyr::pivot_wider(
      names_from = "sample",
      values_from = c("scaled_reads")
    )
  
  # Remove the "lineage_depth" column
  tb <- tb[, !(colnames(tb) %in% c("lineage_depth"))]
  # Reinsert '|' to and move lineages to row names again
  tb <- tb %>% 
    mutate(lineage = str_replace_all(lineage, ">", "\\|"))
  tb <- column_to_rownames(tb, var = "lineage")
}


# This function is called from main.R to perform differential abundance analysis
export("kraken_differential_abundance")
kraken_differential_abundance <- function (kraken_matrix_dir,
                                           metadata_filepath,
                                           da_dir,
                                           statistical_analyses,
                                           css_percentile=0.5) {
  # File input paths
  kraken_analytical <- Sys.glob(glue("{kraken_matrix_dir}/krakenAnalytical_*.csv"))

  metadata <- read.csv(metadata_filepath, header=T)
  
  
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
    map(~ cumNorm(.x, p = css_percentile))
  
  clade_lineages <- row.names(MRcounts(kraken_new_mr[[1]]))
  taxon_lineages <- row.names(MRcounts(kraken_new_mr[[2]]))
  
  # taxa in clade list but not in taxon list are those with 0 taxon counts
  zero_count_lineages <- setdiff(clade_lineages, taxon_lineages)
  # extract rows of those taxa, then set counts to NA
  zero_count_rows <- MRcounts(kraken_new_mr[[1]])[zero_count_lineages,]
  zero_count_rows[TRUE] <- NA
  
  # append the zero-clades to taxon count table, which will be used to
  # recalculate clade counts using scaled taxon counts
  new_clade_counts <- rbind(MRcounts(kraken_css[[2]], norm = TRUE), zero_count_rows)
  new_clade_counts <- calc_clade_counts(new_clade_counts)
  
  exportStats(kraken_css[[2]], file=file.path(da_dir, "CSS_taxa_normalization_stats.tsv"))
  
  
  # Extract the normalized counts into data tables for aggregation
  kraken_norm <- 
    kraken_css %>%
    map(~ data.table(MRcounts(.x, norm=T)))
  
  # Replace clade norm with new experiment with recalculated clade counts
  kraken_norm[[1]] <- data.table(new_clade_counts)
  
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
      ~ kraken_norm$taxonReads[!is.na(eval(as.name(.x))) & eval(as.name(.x)) != "NA" , lapply(.SD, sum), by=.x, .SDcols=!1:10]
    ) %>% 
    set_names(nm=tax_levels)
  
  kraken_taxon_raw_summarised <- 
    tax_levels %>%
    map(
      ~ kraken_raw$taxonReads[!is.na(eval(as.name(.x))) & eval(as.name(.x)) != "NA", lapply(.SD, sum), by=.x, .SDcols=!1:10]
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
  # kraken_species_raw <- kraken_raw[!is.na(Species) & Species != "NA", lapply(.SD, sum), by="Species", .SDcols=!1:8]
  # kraken_species_raw_analytic <- newMRexperiment(counts=kraken_species_raw[, .SD, .SDcols=!"Species"])
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
    analytic <- newMRexperiment(counts=x[, .SD, .SDcols=!"lowest"])
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
    ~ meg_functions$melt_dt(MRcounts(.x), .y) # check if some of this is correct or not? 
  ) # getting warning: binding character and factor vector, coercing into character vector
  
  kraken_taxon_raw_melted <- imap_dfr(
    kraken_taxon_raw_analytic,
    ~ meg_functions$melt_dt(MRcounts(.x), .y)
  ) # getting warning: binding character and factor vector, coercing into character vector
  
  kraken_clade_norm_melted <- imap_dfr(
    kraken_clade_norm_analytic,
    ~ meg_functions$melt_dt(MRcounts(.x), .y)
  )
  
  kraken_clade_raw_melted <- imap_dfr(
    kraken_clade_raw_analytic,
    ~ meg_functions$melt_dt(MRcounts(.x), .y)
  ) # getting warning: binding character and factor vector, coercing into character vector
  
  
  # kraken_taxon_norm_melted$Level_ID <- reorder_tax_ranks(kraken_taxon_norm_melted$Level_ID)
  # kraken_taxon_raw_melted$Level_ID <- reorder_tax_ranks(kraken_taxon_raw_melted$Level_ID)
  # kraken_clade_norm_melted$Level_ID <- reorder_tax_ranks(kraken_clade_norm_melted$Level_ID)
  # kraken_clade_raw_melted$Level_ID <- reorder_tax_ranks(kraken_clade_raw_melted$Level_ID)
  
  # Match metadata ----------------------------------------------------------

  matrix_sampleIDs <- colnames(MRcounts(kraken_clade_norm_analytic[[1]]))
  metadata_sampleIDs <- metadata[, "ID"]
  id_match <- match(matrix_sampleIDs, metadata_sampleIDs)

  if (anyNA(id_match)) {
    stop("Sample IDs in metadata table do not match sample IDs from kraken reports")
  }
  metadata <- data.table(metadata[id_match, ])
  
  # Pair metadata with kraken data
  for( l in 1:length(kraken_taxon_norm_analytic) ) {
    sample_idx <- match(colnames(MRcounts(kraken_taxon_norm_analytic[[l]])), metadata[["ID"]])
    pData(kraken_taxon_norm_analytic[[l]]) <- data.frame(
      metadata[sample_idx, .SD, .SDcols=!"ID"])
    rownames(pData(kraken_taxon_norm_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols="ID"][["ID"]]
    fData(kraken_taxon_norm_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_taxon_norm_analytic[[l]])))
    rownames(fData(kraken_taxon_norm_analytic[[l]])) <- rownames(MRcounts(kraken_taxon_norm_analytic[[l]]))
    pData(kraken_taxon_norm_analytic[[l]]@expSummary$expSummary)$normFactors <- calcNormFactors(kraken_new_mr$taxonReads, p=0.5)
  }
  
  # Use the taxon norm factors because scaled clade counts were calculated from scaled taxon counts
  for( l in 1:length(kraken_clade_norm_analytic) ) {
    sample_idx <- match(colnames(MRcounts(kraken_clade_norm_analytic[[l]])), metadata[["ID"]])
    pData(kraken_clade_norm_analytic[[l]]) <- data.frame(
      metadata[sample_idx, .SD, .SDcols=!"ID"])
    rownames(pData(kraken_clade_norm_analytic[[l]])) <- metadata[sample_idx, .SD, .SDcols="ID"][["ID"]]
    fData(kraken_clade_norm_analytic[[l]]) <- data.frame(Feature=rownames(MRcounts(kraken_clade_norm_analytic[[l]])))
    rownames(fData(kraken_clade_norm_analytic[[l]])) <- rownames(MRcounts(kraken_clade_norm_analytic[[l]]))
    pData(kraken_clade_norm_analytic[[l]]@expSummary$expSummary)$normFactors <- calcNormFactors(kraken_new_mr$taxonReads, p=0.5)
  }
  
  kraken_taxon_names <- names(kraken_taxon_raw_analytic)
  kraken_clade_names <- names(kraken_clade_raw_analytic)
  
  # Apply differential abundance analysis
  for (a in 1:length(statistical_analyses)){
    meg_functions$meg_fitZig(data_list=kraken_taxon_norm_analytic,
                             data_names=kraken_taxon_names,
                             metadata=metadata,
                             zero_mod=model.matrix(~1 + log(libSize(kraken_css$taxonReads))),
                             data_mod=statistical_analyses[[a]]$model_matrix,
                             filter_min_threshold=0.15,
                             contrast_list=statistical_analyses[[a]]$contrasts,
                             random_effect_var=statistical_analyses[[a]]$random_effect,
                             outdir=da_dir,
                             analysis_name=statistical_analyses[[a]]$name,
                             analysis_subset=statistical_analyses[[a]]$subsets,
                             data_type="Microbiome_taxonReads",
                             pval=0.1,
                             top_hits=1000)
  }
  
  for (a in 1:length(statistical_analyses)){
    meg_functions$meg_fitZig(data_list=kraken_clade_norm_analytic,
                             data_names=kraken_clade_names,
                             metadata=metadata,
                             zero_mod=model.matrix(~1 + log(libSize(kraken_css$cladeReads))),
                             data_mod=statistical_analyses[[a]]$model_matrix,
                             filter_min_threshold=0.15,
                             contrast_list=statistical_analyses[[a]]$contrasts,
                             random_effect_var=statistical_analyses[[a]]$random_effect,
                             outdir=da_dir,
                             analysis_name=statistical_analyses[[a]]$name,
                             analysis_subset=statistical_analyses[[a]]$subsets,
                             data_type="Microbiome_cladeReads",
                             pval=0.1,
                             top_hits=1000)
  }
}