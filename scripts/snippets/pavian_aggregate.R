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


widen_results_function <- function(krakenReportPaths, outdir) {
  
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
  
  merged_wide <- krakenReportsPavianMerged %>% 
    select(-percentage) %>% 
    pivot_wider(names_from = Sample, 
                values_from = c(cladeReads, taxonReads),
                values_fill = 0) %>% 
    relocate(name)
  
  # make_kraken_analytical <- function(x, column){
  #   x <- x %>%
  #     select(Sample, name, taxRank, taxID, column, taxLineage) %>%
  #     pivot_wider(names_from = Sample, values_from = column, values_fill = 0) %>%
  #     rename(Lineage = taxLineage)
  #   x
  # }
  # 
  # tax_columns <- c("cladeReads", "taxonReads")
  # 
  # krakenAnalytical <- map(tax_columns, ~ make_kraken_analytical(krakenReportsPavianMerged, .x)) %>%
  #   set_names(nm = tax_columns)
  
  # # Write matrices
  # iwalk(krakenAnalytical,
  #       ~ write.csv(.x, glue("{outdir}/krakenAnalytical_{.y}.csv"),
  #                   row.names = FALSE))
  
  write.csv(merged_wide, glue("{outdir}/krakenAnalytical.csv"), row.names = FALSE)
}

krakenReportPaths <- c("../data/cor_2020/kraken_reports/COR01e_report.txt",
           "../data/cor_2020/kraken_reports/COR01u_report.txt")
outdir <- "results/cor_2020_test" 
widen_results_function(krakenReportPaths, outdir)
