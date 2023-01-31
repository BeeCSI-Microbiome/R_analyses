# For collection of raw taxonomic classification reads for activity 2. Raw read
# count and relative abundance (within sample) tables to be collected for
# sharing with collaborators

packages <- c("tidyverse",
              "vegan",
              "modules",
              "data.table",
              "ggplot2",
              "glue",
              "magrittr")
lapply(packages, library, character.only = TRUE)

dataset_names <- c("cac_2020",
                   "cac_2021",
                   "cas_2020",
                   "cas_2021",
                   "cra_2020",
                   "cra_2021",
                   "hbb_2020",
                   "hbb_2021",
                   "soy_2020",
                   "cor_2020",
                   "lbb_2021",
                   "app_2021")

file_paths <- glue("../data/{dataset_names}/{dataset_names}_aggregated_counts.csv")

# filter out host (apis mellifera), write filtered table
format_tables <- function(file_path) {
  raw_table <- read_csv(file_path)

  # get only timepoint 2 for tables with timepoints
  if (any(str_detect(names(raw_table), "_t2$"))) {
    raw_table <- raw_table %>%
      select(-matches("_t[134]$"))
  }
  # handle clade table
  clade_table <- raw_table %>%
    select(-contains("taxonReads"))

  apis_read_counts <- clade_table %>%
    filter(name == "Apis mellifera") %>%
    select(contains("cladeReads")) %>%
    pivot_longer(cols = everything(), names_to = "sample", values_to = "apis_count")
  apis_tax_lineage <- filter(clade_table, name == "Apis mellifera")$taxLineage

  clade_table2 <- clade_table %>%
    pivot_longer(cols = contains("cladeReads"),
                 names_to = "sample",
                 values_to = "read_count") %>%
    merge(apis_read_counts) %>%
    # if target lineage is a substring of apis_lineage, remove apis count from
    # target taxa
    mutate(read_count = if_else(str_detect(apis_tax_lineage, taxLineage),
                                read_count - apis_count,
                                read_count)) %>%
    select(-apis_count) %>%
    pivot_wider(names_from = sample, values_from = read_count)

  # drop rows containing all 0s
  clade_table2 <- filter(clade_table2, !if_all(contains("cladeReads"), ~ . == 0))
  clade_table2[is.na(clade_table2)] <- 0

  clade_output_path <- basename(file_path) %>%
    str_extract("^\\w{3}_\\d{4}") %>%
    paste0("results/activity2_counts/", ., "_non-host_clade_counts.csv")

  write_csv(clade_table2, clade_output_path)


  # handle taxon counts
  taxon_table <- raw_table %>%
    select(-contains("cladeReads"))

  # Merge with clade table ensuring dropped rows are consistent between tables
  # Drop rows with all NAs (i.e. clade-count-only rows)
  taxon_table2 <- merge(clade_table2, taxon_table) %>%
    select(-contains("cladeReads")) %>%
    filter(!if_all(contains("taxonReads"), is.na))
  taxon_table2[is.na(taxon_table2)] <- 0

  taxon_output_path <- basename(file_path) %>%
    str_extract("^\\w{3}_\\d{4}") %>%
    paste0("results/activity2_counts/", ., "_non-host_taxon_counts.csv")

  write_csv(taxon_table2, taxon_output_path)

  list(clade = clade_table2,
       taxon = taxon_table2)
}

# calculate relative abundance of classified taxa ("unclassified" removed)
make_relative_tables <- function(tb, file_path) {
  taxon_table <- tb$taxon %>%
    filter(!str_detect(name, "^unclassified$|cellular organisms"))

  taxon_table_long <- taxon_table %>%
    pivot_longer(cols = contains("taxonReads"),
                 names_to = "sample",
                 values_to = "read_count")

  total_reads_by_sample <- taxon_table_long %>% group_by(sample) %>%
    summarize(total_reads = sum(read_count))

  relabund_table <- merge(taxon_table_long, total_reads_by_sample) %>%
    mutate(percent = read_count/total_reads*100) %>%
    select(-c(read_count, total_reads)) %>%
    pivot_wider(names_from = sample, values_from = percent)

  relabund_output_path <- basename(file_path) %>%
    str_extract("^\\w{3}_\\d{4}") %>%
    paste0("results/activity2_counts/", ., "_non-host_relative_abundance.csv")

  write_csv(relabund_table, relabund_output_path)

  taxon_table2
}

formatted_tables <- lapply(file_paths, format_tables)

relative_tables <- mapply(make_relative_tables, formatted_tables, file_paths, SIMPLIFY = FALSE)
