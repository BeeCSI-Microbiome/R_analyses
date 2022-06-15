# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

# Identifying "indicator taxa" - those found more often in one treatment group
# compared to another

library(indicspecies)

# Create output dir if it doesn't exist yet
ind_sp_dir <-
  glue("results/{dataset_name}/indicator_species_analysis")
ifelse(!dir.exists(ind_sp_dir),
       dir.create(ind_sp_dir, mode = "777"),
       FALSE)

# calculates proportions and convert to data frame
calc_prop <- function(d) {
  sample_col <- select(d, c("sample", "replicate", "treatment"))
  prop_data <- select(d,-c("sample", "replicate", "treatment")) %>%
    apply(
      MARGIN = 1,
      FUN = function(x)
        x / sum(x) * 100
    ) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample, replicate, treatment)
  
  return(prop_data)
}

format_table <- function(tb, rank_symbol = "all", act_num) {
  sample_regex <- if (act_num == 1) {
    "(.{3,5}(\\d\\d)_d(.)).*"
  } else if (act_num == 2) {
    "(^[[:upper:]]{3}(\\d{2})([u,e])).*"
  } else {
    stop("Activity number must be specified as 1 or 2")
  }
  
  tt <- if (rank_symbol == "all") {
    tb
  } else {
    filter(tb, taxRank == rank_symbol)
  }
  tt <- tt %>%
    pivot_longer(
      cols = contains("Reads"),
      names_to = c("sample", "replicate", "treatment"),
      names_pattern = sample_regex,
      values_to = "read_count"
    ) %>%
    mutate(
      replicate = paste0("Rep ", replicate),
      treatment = paste0("Treatment ", treatment),
      name = case_when(
        name == "Actinobacteria" & taxRank == "C" ~ "Actinobacteria_class",
        TRUE ~ name
      )
    )
  
  tt <- tt %>%
    select(-c(taxRank, lineage))  %>%
    pivot_wider(names_from = "name", values_from = "read_count")
  
  tt <- calc_prop(tt)
}

tt <- format_table(tb = CTX_tables,
                   rank_symbol = "all",
                   act_num = 1)

CAC_tt <- format_table(tb = CAC_tables,
                       rank_symbol = "all",
                       act_num = 2)
