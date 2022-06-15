# Author(s): Lance Lansing
# Produce a summary of read information

# ----------------------------------- Setup ------------------------------------
import("glue")
import("dplyr")
import("stringr")
import("utils")
# ------------------------------------------------------------------------------

# ---------------------------------- Globals -----------------------------------
taxa_list <- c("root", # total classified reads in a sample
               "unclassified",
               "Apis mellifera", # host
               "Bacteria",
               "Nosema apis", 
               "Nosema ceranae",
               "Crithidia mellificae",
               "Lotmaria passim")
# ------------------------------------------------------------------------------


# ------------------------------ Table Formatting ------------------------------

export("create_summary_table")
create_summary_table <- function(read_table, dataset_name){
  
  # Pivot samples to a column, keep only taxa in taxa_list
  ta <- subset(read_table, name %in% taxa_list) %>% 
    select_if(str_detect(names(.), "name|cladeReads")) %>% 
    tidyr::pivot_longer(
      cols = contains("cladeReads"),
      names_to = "sample",
      values_to = "reads"
    )
  
  # sample name formatting
  ta$sample <- gsub("(.*)\\.cladeReads", "\\1", ta$sample)
  
  # Pivot to get taxa as columns
  ta <- ta %>% 
    tidyr::pivot_wider(
      names_from = "name",
      values_from = c("reads")
    )
  
  # some column name formatting
  ta <- ta %>% dplyr::rename(classified = root) %>%
    dplyr::rename_with(~ gsub(" ", "_", .x)) %>%
    dplyr::rename_with(~ str_c(.x, "_reads"), !matches("sample"))
  
  # create total read count column
  ta <- ta %>% 
    dplyr::mutate(ta, total_reads = unclassified_reads+classified_reads,
                  .after = sample) %>%
    # percent classified
    dplyr::mutate(percent_classified = classified_reads/total_reads * 100,
                  .after = classified_reads) %>%
    # percent unclassified
    dplyr::mutate(percent_unclassified = unclassified_reads/total_reads * 100,
                  .after = unclassified_reads)
  
  # create "percent of classified reads" columns for taxa
  ta <- ta %>% dplyr::mutate(across(
    !matches("sample|total_reads|*classified*"),
    ~ .x/classified_reads*100,
    .names = "percent_classified_{gsub('_reads', '', .col)}"
  ))
  
  # Drop read # columns for taxa
  ta <- ta %>% select(!matches("(?<!total|classified|Bacteria|Apis_mellifera)_reads", perl = TRUE))
  ta <- ta %>% relocate("percent_classified_Bacteria", .after = "percent_classified_Apis_mellifera")
  
  write.table(ta, file = glue("results/{dataset_name}/{dataset_name}_read_summary.csv"), quote = F, sep = ",", row.names = F)
  return(ta)
}