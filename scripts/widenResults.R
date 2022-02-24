# This script reads in kraken reports and agglomerates them into analytic 
# matrices for both clade and taxon reads via Pavian.
# Load packages -----------------------------------------------------------

import("pavian")
import("glue")
import("dplyr")
import("purrr")
import("stringr")
import("tidyr")
import("utils")


export("widen_results_function")
widen_results_function <- function(krakenReportPaths, krakenReportNames, outdir) {
  
  krakenReportNames <-
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
        mutate(taxLineage = str_replace(taxLineage, "-_root\\|-_cellular organisms\\|", "")) %>%
        mutate(taxLineage = str_replace(taxLineage, "-_root\\|", ""))
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
    set_names(nm = tax_columns)
  
  # Write matrices
  iwalk(krakenAnalytical,
        ~ write.csv(.x, glue("{outdir}/krakenAnalytical_{.y}.csv"),
                    row.names = FALSE))
}