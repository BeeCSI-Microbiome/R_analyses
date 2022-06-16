# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

# Identifying "indicator taxa" - those found more often in one treatment group
# compared to another

library(indicspecies)

# Create output dir if it doesn't exist yet
ind_sp_dir <- glue("results/{dataset_name}/indicator_species_analysis")
ifelse(!dir.exists(ind_sp_dir),
       dir.create(ind_sp_dir, mode = "777"), FALSE)

# calculates proportions and convert to data frame
calc_prop <- function(d) {
  sample_col <- select(d, c("sample", "replicate", "treatment"))
  prop_data <- select(d, -c("sample", "replicate", "treatment")) %>%
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample, replicate, treatment)
  
  return(prop_data)
}
# pass in filepath to tables file
run_indicator_analysis <- function(table_string, analysis_func, rank_symbol="all") {
  # Pivot tables into wide format (1 column per taxa, 1 row per sample) and
  # create metadata category columns
  if(rank_symbol == "all"){
    tt <- tables[[table_string]]
  } else{
    tt <- filter(tables[[table_string]], taxRank == rank_symbol)
  }
  tt <- tt %>%
    pivot_longer(
      cols = contains("Reads"), # pivot columns containing "Reads"
      names_to = c("sample", "replicate", "treatment"), # make 3 new columns named these. These will contain names of tt cols
      names_pattern = "(CTX(\\d*)_d(.)).*", # split up tt columns names to bin parts in names_to
      values_to = "read_count"
    ) %>%
    mutate(
      treatment = case_when(
        treatment == "0" ~ "Control",
        treatment == "C" ~ "Clothianidin",
        treatment == "T" ~ "Thiamethoxam"
      ),
      replicate = paste0("Rep ", replicate),
      name = case_when(
        name == "Actinobacteria" & taxRank == "C" ~ "Actinobacteria_class",
        TRUE ~ name
      )
    )
  # tt <- if(rank_symbol == "all") { tt }
  # else { filter(tt, taxRank == rank_symbol) }
  # browser()
  tt <- tt %>% 
    select(-c(taxRank, lineage))  %>%
    pivot_wider(names_from = "name", values_from = "read_count")
  
  
  tt <- calc_prop(tt)
  
  # Extract count table and metadata vectors
  tt_abund <- tt[,4:ncol(tt)] # generalize maybe
  tt_treat <- tt$treatment
  tt_rep <- tt$replicate
  
  # Analysis function string for output naming
  func_out_str <- str_replace(analysis_func, "\\.", "")
  
  # Analyze treatment and replicates. Use sink() for writing summaries to files
  inv_treat = multipatt(tt_abund, tt_treat, func = analysis_func, control = how(nperm=9999))
  inv_rep = multipatt(tt_abund, tt_rep, func = analysis_func, control = how(nperm=9999))
  sink(glue("{ind_sp_dir}/{func_out_str}_{table_string}_treatment_indicators_{rank_symbol}.txt"))
  summary(inv_treat)
  sink()
  sink(glue("{ind_sp_dir}/{func_out_str}_{table_string}_replicate_indicators_{rank_symbol}.txt"))
  summary(inv_rep)
  sink()
}


run_indicator_analysis("raw_taxon", "r.g")
run_indicator_analysis("raw_clade", "r.g")

run_indicator_analysis("raw_taxon", "r.g", "S")
# run_indicator_analysis("scaled_taxon", "r.g", "S")
run_indicator_analysis("raw_clade", "r.g", "S")
# run_indicator_analysis("scaled_clade", "r.g", "S")

run_indicator_analysis("raw_taxon", "r.g", "G")
# run_indicator_analysis("scaled_taxon", "r.g", "G")
run_indicator_analysis("raw_clade", "r.g", "G")
# run_indicator_analysis("scaled_clade", "r.g", "G")
