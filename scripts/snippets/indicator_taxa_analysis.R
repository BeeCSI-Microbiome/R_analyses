######################################################################
## Indicator Taxa Analysis
######################################################################
## Identifying "indicator taxa" - those found more 
## often in one treatment group compared to another
##
## Using:
## https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html
######################################################################


setwd("~/Github/R_analyses")

library(indicspecies)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(glue)




#####################################################
## takes in activity type
## returns the regular expression to use in pivoting
#####################################################
reg_ex <- function(activity){
  if(activity == 1){
    reg <- paste0("((?:(?:\\w_){0,1}\\w{3})(\\d\\d)_d(\\d)(?:_t(\\d)){0,1}).*")
  } 
  if(activity == 2){
    reg <- paste0("(\\D\\D\\D(..)([e|u])(?:_t(\\d)){0,1}).*")
  }
  return(reg)
}

#####################################################
## takes in the regular expression being used and one
##    sample id from a column in tt
## returns the metadata columns to use when pivoting
#####################################################
col_names <- function(reg_ex, sampleID){
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

#####################################################
## takes in the tt table from run indicator analysis
##    and the metadata columns being used
## returns the proportions for each taxa as a table
#####################################################
calc_prop <- function(tt, meta_cols) {
  sample_col <- select(tt, all_of(meta_cols))
  prop_data <- as.matrix(select(tt, -all_of(meta_cols))) %>%
    apply(MARGIN = 1, FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(all_of(meta_cols))
  
  return(prop_data)
}

#####################################################
## takes in data set name (same as in main.R), the 
##    activity type of the data set, and the string 
##    for the desired table (e.g. "raw_taxon"), the 
##    analysis function you'd like to use, and the
##    rank_symbol you'd like
## Does the indicator taxa analysis and writes those 
##    results along with adjusted p-values to 2 files
##    (one for treatment and one for replicate)
#####################################################
run_indicator_analysis <- function(dataset_name, activity, table_string, analysis_func, rank_symbol="all") {
  # Create output dir if it doesn't exist yet
  ind_sp_dir <- glue("results/{dataset_name}/indicator_species_analysis")
  ifelse(!dir.exists(ind_sp_dir),
         dir.create(ind_sp_dir, mode = "777"), FALSE)
  
  # read in tables files from the relevant results folder
  file <- paste0("results/", dataset_name,"/", table_string,"_counts.csv")
  table <- read.csv(file)
  if(rank_symbol == "all"){
    tt <- table
  } else{
    tt <- filter(table, taxRank == rank_symbol)
  }
  
  # Pivot tables into wide format (1 column per taxa, 1 row per sample) and
  # create metadata category columns
  d_nm <- substr(dataset_name, 1, nchar(dataset_name)-5)
  d_ctl <- paste0(substring(d_nm, 1, 2), "ctl")
  regx <- reg_ex(activity)
  sample_id <- colnames((tt %>% select(-c("name", "taxRank", "taxID", "depth", "taxLineage")))[1])
  meta_cols <- col_names(regx, sample_id)
  tt <- tt %>%
    pivot_longer(
      cols = contains(c(d_nm, d_ctl)), 
      names_to = meta_cols,
      names_pattern = regx, 
      values_to = "read_count"
    )
  
  # make metadata columns more readable
  tt$treatment[-which(tt$treatment=="0")] <-  paste0("Treatment ", tt$treatment[-which(tt$treatment=="0")])
  tt$treatment[which(tt$treatment=="e")] <- "Treatment"
  tt$treatment[which(tt$treatment=="0" | tt$treatment=="u")] <- "Control"
  tt$replicate <- paste0("Rep ", tt$replicate)
  
  # Add taxa IDs to the names in order to deal with duplicate names 
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
  
  # Analysis function string for output naming
  func_out_str <- str_replace(analysis_func, "\\.", "")
  
  # Analyze treatment and replicates
  inv_treat = multipatt(tt_abund, tt_treat, func = analysis_func, control = how(nperm=9999))
  inv_rep = multipatt(tt_abund, tt_rep, func = analysis_func, control = how(nperm=9999))
  
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
  sink(glue("{ind_sp_dir}/{func_out_str}_{table_string}_treatment_indicators_{rank_symbol}.txt"))
  summary(inv_treat)
  cat("\n\n\nAdjusted p-values\n")
  cat("-----------------\n")
  print(inv_treat.sign)
  sink()
  sink(glue("{ind_sp_dir}/{func_out_str}_{table_string}_replicate_indicators_{rank_symbol}.txt"))
  summary(inv_rep)
  cat("\n\n\nAdjusted p-values\n")
  cat("-----------------\n")
  print(inv_rep.sign)
  sink()
}







# dataset to run analysis on
dataset_name <- "pdv_2021"
activity <- 1

# run for all taxa
run_indicator_analysis(dataset_name, activity, "raw_taxon", "r.g")
run_indicator_analysis(dataset_name, activity, "raw_clade", "r.g")
# run for species taxa
run_indicator_analysis(dataset_name, activity, "raw_taxon", "r.g", "S")
run_indicator_analysis(dataset_name, activity, "raw_clade", "r.g", "S")
# run for genus taxa
run_indicator_analysis(dataset_name, activity, "raw_taxon", "r.g", "G")
run_indicator_analysis(dataset_name, activity, "raw_clade", "r.g", "G")

