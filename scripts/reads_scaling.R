# Author(s): Lance Lansing
# Oct 2021
# Normalization via Cumulative Sum Scaling of read count data
# Input: table from Pavian with clades AND taxon counts, not collapsed


# ---------------------------------- Imports -----------------------------------
import('ggplot2')
import('stringr')
import('metagenomeSeq')
import('dplyr')
import('stats', 'setNames', 'aggregate')
# ______________________________________________________________________________


export("format_count_table")
format_count_table <- function(tb){
  # Performs some preliminary formatting on count table
  
  # Remove the "Max" columns
  tb <- tb[, !(colnames(tb) %in% c("Max", "Max.1"))]
  
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


export('scaling_procedure')
scaling_procedure <- function(tb, css_percentile){
  # Given the formatted and filtered table, return a list of four ordered items:
  # (Raw taxon table, raw clade table, scaled taxon table, scaled clade table)
  
  tb_raw_taxon <- drop_all_NA_rows(tb)
  # Get clade counts. We calculate counts rather than use clade counts from Pavian
  # in order to account for taxa that we filtered out
  tb_raw_clade <- get_clade_data(tb)
  
  # Perform the scaling
  tb_scaled_list <- css_scale(tb, css_percentile)
  
  # Get scaled taxon table by dropping taxa rows with no counts in any sample
  tb_scaled_taxon <- drop_all_NA_rows(tb_scaled_list$tb_scaled)
  
  # Calculate scaled clade counts from scaled taxon counts 
  tb_scaled_clade <- get_clade_data(tb_scaled_list$tb_scaled)
  
  list(css_MRexp=tb_scaled_list$css_MRexp,
       raw_taxon=tb_raw_taxon,
       raw_clade=tb_raw_clade,
       scaled_taxon=tb_scaled_taxon,
       scaled_clade=tb_scaled_clade)
}


export('drop_all_NA_rows')
drop_all_NA_rows <- function(tb){
  # drops all rows from the table in which counts are NA for all samples
  counts_only <- as.data.frame(select(tb, contains('taxonReads')))
  tb <- tb[as.logical((rowSums(is.na(counts_only))-ncol(counts_only))), ]
  tb[is.na(tb)] <- 0
  tb
}


export('get_clade_data')
get_clade_data <- function(tb){
  # We calculate counts rather than use clade counts from Pavian in order to 
  # account for taxa that we filtered out
  tb <- calc_clade_counts(tb)
  names(tb) <- gsub('(.*).taxonReads', '\\1.cladeReads', names(tb))
  tb
}


export('css_scale')  
css_scale <- function(tb, css_percentile){

  # get table with only counts
  counts_only <- as.data.frame(select(tb, ends_with('taxonReads')))
  
  # Get table of raw counts with only the taxa that are not ALL NA
  # (ie filter the taxa with clade counts but no taxon counts)
  # This table will be used in sum visualization later
  taxon_raw <- tb[as.logical((rowSums(is.na(counts_only))-ncol(counts_only))), ]
  
  # set row names
  row.names(counts_only) <- tb$lineage
  # remove rows with all NA
  counts_only <- counts_only[as.logical((rowSums(is.na(counts_only))-ncol(counts_only))), ]
  # Set remaining NA values to 0
  counts_only[is.na(counts_only)] <- 0
  
  css_MRexp <- newMRexperiment(counts_only)
  
  # Perform CSS normalization
  css_MRexp <- cumNorm(css_MRexp, p = css_percentile)
  
  # Get scaled counts into table
  scaled_counts <- data.frame(MRcounts(css_MRexp, norm = TRUE))
  scaled_counts <- tibble::rownames_to_column(scaled_counts, var = 'lineage')
  
  # Merge scaled_counts back with tb, overwriting count values
  tb <- merge(tb, scaled_counts, by = 'lineage', all = T)
  # Remove the unnormalized columns
  tb <- select(tb, !ends_with('taxonReads.x'))
  names(tb) <- gsub(names(tb), pattern='(.*).y', replacement='\\1')
  
  # Another function, writes a visualization  'results/scaling_visualization.png'
  visualize_scaling(taxon_raw, scaled_counts, css_percentile)
  
  # return table
  list(tb_scaled=tb, css_MRexp=css_MRexp)
}


visualize_scaling <- function(taxon_raw, scaled_counts, css_percentile){
  # Prepare normalized data -
  taxon_norm <- select(merge(taxon_raw, scaled_counts, by = 'lineage', all = T),
                       !ends_with('taxonReads.x'))
  
  # Pivot samples into 1 column
  taxon_norm <- taxon_norm %>% 
    tidyr::pivot_longer(
      cols = contains('taxonReads'),
      names_to = "sample",
      values_to = "scaled_reads"
    )
  
  # Clean up the sample name entries
  taxon_norm$sample <- gsub('(\\w+\\d\\d(u|e)).*$', "\\1", taxon_norm$sample)
  
  # Aggregate to get read sums by sample
  taxon_norm <- setNames(aggregate(taxon_norm$scaled_reads,
                                   by = list(sample = taxon_norm$sample),
                                   FUN = sum),
                         c('sample','scaled_reads'))
  
  # Prepare raw data
  # Pivot samples into 1 column
  taxon_raw <- taxon_raw %>% 
    tidyr::pivot_longer(
      cols = contains('taxonReads'),
      names_to = "sample",
      values_to = "reads"
    )
  
  # Clean up the sample name entries
  taxon_raw$sample <- gsub('(\\w+\\d\\d(u|e)).*$', "\\1", taxon_raw$sample)
  # replace NA with 0s
  taxon_raw[is.na(taxon_raw)] <- 0
  # Aggregate to get sums by sample
  taxon_raw <- setNames(aggregate(taxon_raw$reads,
                                  by = list(sample = taxon_raw$sample), FUN = sum),
                        c('sample','reads'))
  
  # Merge raw/scaled
  taxon_all <- merge(taxon_norm, taxon_raw)
  
  # Get proportions by sample
  read_sum <- sum(taxon_all$reads)
  norm_read_sum <- sum(taxon_all$scaled_reads)
  
  # Pivot raw and scaled into 1 column
  taxon_all <- taxon_all %>% 
    tidyr::pivot_longer(
      cols = 2:3,
      names_to = "scaled",
      values_to = "reads"
    )
  # entry formatting
  taxon_all$scaled <- ifelse(str_detect(taxon_all$scaled, "scaled"), T, F)
  
  # Create proportions by sample using the respective sums for scaled vs raw
  taxon_all$proportion <- ifelse(taxon_all$scaled,
                                 (taxon_all$reads / norm_read_sum),
                                 (taxon_all$reads / read_sum))
  
  # Visualize with bar plot
  # TODO: Titles may have to change, depending on scaling decisions 
  scaling_vis <- taxon_all %>%
    ggplot(aes(x = sample, y = proportion, fill = scaled)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Total read proportion by sample, before and after cumulative sum scaling",
         subtitle = sprintf("All taxon reads, A. mellifera and unclassified reads filtered, CSS percentile = %s", as.character(css_percentile))) +
    geom_text(aes(label = round(proportion*100, 1)), vjust = 1.6, color = "white",
              position = position_dodge(0.9), size = 3.5)+
    scale_fill_brewer(palette = "Paired")+
    theme_minimal()
  
  scaling_vis
  
  ggsave(plot = scaling_vis, filename = 'results/scaling_visualization.png', bg = 'white')
}

export('calc_clade_counts')
calc_clade_counts <- function(tb){
  # Reaggregate clade counts with normalized values
  
  # Count lineage_depths. Clade counts will be counted from leaves up
  tb <- tb %>%
    mutate(lineage_depth = str_count(lineage, '>'))
  
  # Set remaining NA values to 0
  tb[is.na(tb)] <- 0
  
  # Pivot samples into 1 column
  tb <- tb %>% 
    tidyr::pivot_longer(
      cols = contains('taxonReads'),
      names_to = "sample",
      values_to = "scaled_reads"
    )
  
  # Clean sample names
  tb$sample <- gsub('(.*)\\.taxonReads.y', '\\1', tb$sample)
  
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
      filter(tb,
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
  tb <- tb %>% relocate(lineage, .after = everything())
}
