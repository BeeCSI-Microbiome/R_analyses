# Author(s): Lance Lansing
# Nov 2021
#  Visualize the Cumulative Sum Scaling
#
# This file is called with source() from ANOTHER FILE. As such, it expects
# certain objects and data tables to be in the R workspace


# Prepare normalized data ----------------------------------------------------
taxon_norm <- select(merge(taxon_raw, normd, by = 'lineage', all = T), !ends_with('taxonReads.x'))

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
                                 by = list(sample = taxon_norm$sample), FUN = sum),
                       c('sample','scaled_reads'))


# Prepare raw data -----------------------------------------------------------
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


# Merge raw/scaled -----------------------------------------------------------
taxon_all <- merge(taxon_norm, taxon_raw)


# Get proportions by sample --------------------------------------------------
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


# Visualize with bar plot ----------------------------------------------------
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
  # ______________________________________________________________________________