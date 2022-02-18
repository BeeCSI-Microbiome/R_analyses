# This script reads in kraken reports and agglomerates them into an analytic 
# matrix with both clade and taxon reads via Pavian.
# Load packages -----------------------------------------------------------

import("tidyverse")
import("pavian")
import("here")

# Create directory for aggregated data ------------------------------------

#aggregated_dir <- here("aggregated_data_for_analysis")

#ifelse(
 # !dir.exists(aggregated_dir), 
 # dir.create((aggregated_dir), mode='777'), 
 # FALSE
#)

# Read Kraken reports with Pavian ---------------------------------------------

# Obtain the filenames of all Kraken reports

#krakenReportPaths <- Sys.glob("../data/hbb_2020/kraken_reports/*_report.txt")

#krakenReportNames <- list.files(path = "../data/_2020/kraken_reports/", pattern = '*_report.txt')
export("widen_results_function")

widen_results_function <- function(krakenReportPaths) {

krakenReportNames <- #does this need to be called? 
  krakenReportNames %>%
  map(function(x) str_replace(x, "\\_report.txt$", ""))

# Read in reports via Pavian
krakenReportsPavian <-
  krakenReportPaths %>%
  map(function(x) pavian::read_report(x)) %>%
  set_names(nm = krakenReportNames)

krakenReportsPavian <-
  krakenReportsPavian %>%
  map(safely(function(x){
    filt <- pavian::filter_taxon(
      report = x, 
      filter_taxon = "Eukaryota", 
      rm_clade = TRUE, 
      do_message = TRUE
    )
    filt
  }))

# Remove some top-level meta groupings
taxa_to_remove <- c("u_unclassified", "-_root", "-_cellular organisms")

# Merge reports
krakenReportsPavianMerged <- 
  krakenReportsPavian %>%
  map_dfr(function(x){
    x <- x$result %>%
      filter(!name %in% taxa_to_remove) %>%
      filter(taxRank != "-") %>%
      mutate(taxLineage=str_replace(taxLineage, "-_root\\|-_cellular organisms\\|", "")) %>%
      mutate(taxLineage=str_replace(taxLineage, "-_root\\|", ""))
    x
  }, .id = "Sample")

make_kraken_analytical <- function(x, column){
  x <- x %>%
    select(Sample, column, taxLineage) %>%
    spread(key = Sample, value = column, fill = 0) %>%
    rename(Lineage = taxLineage)
  x
}

tax_columns <- c("cladeReads", "taxonReads")

krakenAnalytical <- map(tax_columns, ~ make_kraken_analytical(krakenReportsPavianMerged, .x)) %>%
  set_names(nm=tax_columns)

# Write matrices
iwalk(krakenAnalytical,
      ~ write.csv(.x, here("aggregated_data_for_analysis", paste0("krakenAnalytical", "_",.y,".csv")), row.names = FALSE))
}