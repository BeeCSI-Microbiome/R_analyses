# Author(s): Lance Lansing and Emma Lee
# Summer 2022
# Identifying "indicator taxa" - those found more 
# often in one treatment group compared to another
# 
# Using:
# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html


# ---------------------------------- Imports -----------------------------------
import('indicspecies')
import('data.table')
import('dplyr')
import('tidyr')
import('stringr')
import('glue')
import('permute')
import('stats')
# ______________________________________________________________________________

taxa_lvl_key <- c(D="Domain",
                  P="Phylum",
                  C="Class",
                  O="Order",
                  F="Family",
                  G="Genus",
                  S="Species",
                  all="all")

# --------------------------------- Functions ----------------------------------
export("run_indicator_analysis")
run_indicator_analysis <- function(dataset_name, table, ind_sp_dir, treatment_key, rank_symbol="all") {
  ## takes in data set name (same as in main.R), the 
  ##    string for the desired table (e.g. "raw_taxon"),
  ##    the analysis function you'd like to use, and 
  ##    the rank_symbol you'd like
  ## Does the indicator taxa analysis and writes those 
  ##    results along with adjusted p-values to 2 files
  ##    (one for treatment and one for replicate)
  
  
  if(rank_symbol == "all"){
    tt <- table
  } else{
    tt <- dplyr::filter(table, taxRank == rank_symbol)
  }
  
  # Pivot tables into wide format (1 column per taxa, 1 row per sample) and
  # create metadata category columns
  sample_id <- colnames((tt %>% select(-c("name", "taxRank", "taxID", "depth", "taxLineage")))[1])
  regx <- reg_ex(sample_id)
  meta_cols <- col_names(regx, sample_id)
  tt <- tt %>%
    pivot_longer(
      cols = !c("name", "taxRank", "taxID", "depth", "taxLineage"), 
      names_to = c("sample", "replicate", "treatment", "time"),
      names_pattern = regx, 
      values_to = "read_count"
    )
  if(!("time" %in% meta_cols)){
    tt <- tt %>% select(-c(time))
  }
  
  # make metadata columns more readable
  tt$treatment <- sapply(c(1:length(tt$treatment)), function(x){
    treatment_key[[which(names(treatment_key)==tt$treatment[x])]]
  })
  tt$treatment[which(tt$treatment=="e")] <- "Treatment"
  tt$treatment[which(tt$treatment=="u")] <- "Control"
  tt$replicate <- paste0("Rep ", tt$replicate)

  name_and_id <- paste0(tt$name, " (", tt$taxID, ")")
  tt <- cbind(name_and_id, tt)
  
  tt <- tt %>% 
    select(-c(taxRank, taxLineage, taxID, depth, name))  %>%
    pivot_wider(names_from = "name_and_id", values_from = "read_count") 
  
  # calculate proportions
  tt <- calc_prop(tt, meta_cols)
  
  # Extract count table and metadata vectors
  tt_abund <- tt[,(length(meta_cols)+1):ncol(tt)] 
  tt_treat <- tt$treatment
  tt_rep <- tt$replicate

  
  # Analyze treatment and replicates
  inv_treat = multipatt(tt_abund, tt_treat, func = "r.g" , control = how(nperm=9999))
  inv_rep = multipatt(tt_abund, tt_rep, func = "r.g", control = how(nperm=9999))
  
  # adjusted p-vals
  #extract table of stats
  inv_treat.sign <- as.data.table(inv_treat$sign, keep.rownames = TRUE)
  #add adjusted p-value
  inv_treat.sign[, p.value.bh := p.adjust(p.value, method = "BH")]
  #now can select only the indicators with adjusted significant p-values
  inv_treat.sign[p.value.bh <= 0.05, ]
  #for replicates now
  inv_rep.sign <- as.data.table(inv_rep$sign, keep.rownames = TRUE)
  #add adjusted p-value
  inv_rep.sign[, p.value.bh := p.adjust(p.value, method = "BH")]
  #now can select only the indicators with adjusted significant p-values
  inv_rep.sign[p.value.bh <= 0.05, ]
  
  # write to file using sink()
  sink(glue("{ind_sp_dir}/{dataset_name}_treatment_indicators_{taxa_lvl_key[rank_symbol]}.txt"))
  summary(inv_treat)
  cat("\n\n\nAdjusted p-values\n")
  cat("-----------------\n")
  print(inv_treat.sign)
  sink()
  sink(glue("{ind_sp_dir}/{dataset_name}_replicate_indicators_{taxa_lvl_key[rank_symbol]}.txt"))
  summary(inv_rep)
  cat("\n\n\nAdjusted p-values\n")
  cat("-----------------\n")
  print(inv_rep.sign)
  sink()
}


reg_ex <- function(sampleID){
  ## takes in sample id and figures out activity 
  ##    number from it
  ## returns the regular expression to use in pivoting
  
  
  if (all(str_detect(sampleID, "_d[[:alnum:]]+"))) {
    # Does sample string match activity 1 pattern?
    reg <- paste0("((?:(?:\\w_){0,1}\\w{3})(\\d\\d)_(d\\d)(?:_t(\\d)){0,1}).*")
  } else if (all(str_detect(sampleID, "\\D\\D\\D[[:alnum:]]{2}[e|u]"))) {
    # Or activity 2 pattern?
    reg <- paste0("(\\D\\D\\D(..)([e|u])(?:_t(\\d)){0,1}).*")
  } else {
    stop("The sample column names do not match the implemented regex patterns")
  }
  return(reg)
}


col_names <- function(reg_ex, sampleID){
  ## takes in the regular expression being used and one
  ##    sample id from a column in tt
  ## returns the metadata columns to use when pivoting
  
  
  matches <- str_match(sampleID, reg_ex)
  matches <- matches[!is.na(matches)]
  if(length(matches) == 5){
    cols <- c("sample", "replicate", "treatment", "time")
  }
  if(length(matches) == 4){
    cols <- c("sample", "replicate", "treatment")
  }
  return(cols)
}


calc_prop <- function(tt, meta_cols) {
  ## takes in the tt table from run indicator analysis
  ##    and the metadata columns being used
  ## returns the proportions for each taxa as a table
  
  
  sample_col <- select(tt, all_of(meta_cols))
  prop_data <- as.matrix(select(tt, -all_of(meta_cols))) %>%
    apply(MARGIN = 1, FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(all_of(meta_cols))
  
  return(prop_data)
}

