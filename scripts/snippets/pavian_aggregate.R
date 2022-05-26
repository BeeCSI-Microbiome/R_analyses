# This script reads in kraken reports and agglomerates them into analytic 
# matrices for both clade and taxon reads via Pavian.
# Load packages -----------------------------------------------------------

packages <- c("pavian",
              "dplyr",
              "purrr",
              "stringr",
              "tidyr",
              "glue")
lapply(packages, library, character.only = TRUE)


widen_results_function <- function(krakenReportPaths, outpath) {
  
  krakenReportNames <-
    krakenReportPaths %>%
    map(function(x) str_replace(basename(x), "\\_report.txt$", ""))
  
  # Read in reports via Pavian
  krakenReportsPavian <-
    krakenReportPaths %>%
    map(function(x) pavian::read_report(x)) %>%
    set_names(nm = krakenReportNames)
  
  # Merge reports
  krakenReportsPavianMerged <- 
    krakenReportsPavian %>%
    map_dfr(function(x){
      x
    }, .id = "Sample")
  
  # Fix some taxa name and lineage formatting
  krakenReportsPavianMerged <- krakenReportsPavianMerged %>% 
    mutate(name = str_replace(name, "^(\\w|-)_", ""),
           taxLineage = str_replace_all(taxLineage, "\\|", ">"),
           taxLineage = str_replace_all(taxLineage, "._", ""))
  
  krakenReportsPavianMerged$taxonReads[krakenReportsPavianMerged$taxonReads == 0] <- NA
  
  merged_wide <- krakenReportsPavianMerged %>% 
    select(-percentage) %>% 
    pivot_wider(names_from = Sample, 
                values_from = c(cladeReads, taxonReads)) %>% 
    relocate(name)
  
  write.csv(merged_wide, outpath, row.names = FALSE)
}

krakenReportPaths <- Sys.glob("../data/ctx_2020/kraken_reports/*d0*_report.txt") %>% 
  c(Sys.glob("../data/ctx_2020/kraken_reports/*dC*_report.txt"))

outpath <- "../data/ctx_2020/CTX_dC_all_taxa.csv" 
widen_results_function(krakenReportPaths, outpath)
