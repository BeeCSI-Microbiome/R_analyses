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
