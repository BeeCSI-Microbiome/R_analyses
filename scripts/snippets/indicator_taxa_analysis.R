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





####################################################
## takes in the tt table from run_indicator_analysis
## returns the type of sample id present 
####################################################
id_type <- function(tt){
  
  # 1) ABC{*}{e,u}_t{*}
  # 2) ABC{*}_d{*}
  # 3) A_CTL{*}_d0 and A_BCD{*}_d{*}
  # 4) ABC{*}_d{*}_t{*}
  # 5) ABC{*}{e,u}
  
  nt <- tt %>%
    select(-c(name, taxRank, taxID, depth, taxLineage))
  tester <- colnames(nt)[1]
  parts <- unlist(str_split(tester, ""), use.names=FALSE)
  counts <- as.data.frame(table(parts[1:3]))
  if("_" %in% counts$Var1){
    which1 <- "3"
  } else{
    counts2 <- as.data.frame(table(parts))
    how_many2 <- filter(counts2, parts == "_")[,2]
    if(identical(how_many2, integer(0))){
      which1 <- "5"
    } else if(how_many2 == 2){
      which1 <- "4"
    } else{
      num_letters <- str_count(tester, "\\D") - how_many
      if(num_letters == 5){
        which1 <- "1"
      }
      if(num_letters == 4){
        which1 <- "2"
      }
    }
  }
  return(which1)
}

#####################################################
## takes in the tt table from run_indicator_analysis,
##    the dataset name, and the type of sample id
##    present
## returns the regular expression to use in pivoting
#####################################################
reg_ex <- function(d_name, id_num, tt){
  if(id_num == 1){
    reg <- paste0("(\\D\\D\\D(..)([e|u])_t(.)).*")
  } 
  if(id_num == 2){
    reg <- paste0("(\\D\\D\\D(..)_d(.)).*")
  }
  if(id_num == 3){
    reg <- "((?:(?:\\w_){0,1}\\w{3})(\\d\\d)_d(\\d)).*"
  } 
  if(id_num == 4){
    reg <- paste0("(\\D\\D\\D(..)_d(.)_t(.)).*")
  }
  if(id_num == 5){
    reg <- paste0("(\\D\\D\\D(..)([e|u])).*")
  }
  return(reg)
}

#####################################################
## takes in the the type of sample id present
## returns the metadata columns to use when pivoting
#####################################################
col_names <- function(id_num){
  if(id_num == 1 | id_num == 4){
    cols <- c("sample", "replicate", "treatment", "time")
  }
  if(id_num == 2 | id_num == 3 | id_num ==5){
    cols <- c("sample", "replicate", "treatment")
  }
  return(cols)
}

#####################################################
## takes in the tt table from run indicator analysis
##    and the sample id type present
## returns the proportions for each taxa as a table
#####################################################
calc_prop <- function(tt, type) {
  meta_cols <- col_names(type)
  sample_col <- select(tt, all_of(meta_cols))
  prop_data <- as.matrix(select(tt, -meta_cols)) %>%
    apply(MARGIN = 1, FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(all_of(meta_cols))
  
  return(prop_data)
}

#####################################################
## takes in dataset name (same as in main.R), the 
##    string for the desired table (e.g. "raw_taxon"),
##    the analysis function you'd like to use, and
##    the rank_symbol you'd like
## Does the indicator taxa anlysis and writes those 
##    results along with adjusted p-values to 2 files
##    (one for treatment and one for replicate)
#####################################################
run_indicator_analysis <- function(dataset_name, table_string, analysis_func, rank_symbol="all") {
  # Create output dir if it doesn't exist yet
  ind_sp_dir <- glue("results/{dataset_name}/indicator_species_analysis")
  ifelse(!dir.exists(ind_sp_dir),
         dir.create(ind_sp_dir, mode = "777"), FALSE)
  
  # read in tables files from the relevant results folder
  file <- paste0("results/",dataset_name,"/",table_string,"_counts.csv")
  table <- read.csv(file)
  if(rank_symbol == "all"){
    tt <- table
  } else{
    tt <- filter(table, taxRank == rank_symbol)
  }
  
  # Pivot tables into wide format (1 column per taxa, 1 row per sample) and
  # create metadata category columns
  d_nm <- substr(dataset_name, 1, nchar(dataset_name)-5)
  type <- id_type(tt)
  if(type==3){
    d_ctl <- paste0(substring(d_nm, 1, 2), "ctl")
  } else{
    d_ctl <- d_nm
  }
  cols <- col_names(type)
  tt <- tt %>%
    pivot_longer(
      cols = contains(c(d_nm, d_ctl)), 
      names_to = cols,
      names_pattern = reg_ex(dataset_name, type, tt), 
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
  tt <- calc_prop(tt, type)
  
  # Extract count table and metadata vectors
  tt_abund <- tt[,(length(cols)+1):ncol(tt)] 
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

# run for all taxa
run_indicator_analysis(dataset_name, "raw_taxon", "r.g")
run_indicator_analysis(dataset_name, "raw_clade", "r.g")
# run for species taxa
run_indicator_analysis(dataset_name, "raw_taxon", "r.g", "S")
run_indicator_analysis(dataset_name, "raw_clade", "r.g", "S")
# run for genus taxa
run_indicator_analysis(dataset_name, "raw_taxon", "r.g", "G")
run_indicator_analysis(dataset_name, "raw_clade", "r.g", "G")

