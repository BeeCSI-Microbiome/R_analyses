# https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html

# Identifying "indicator taxa" - those found more often in one treatment group
# compared to another

setwd("~/Github/R_analyses")
dataset_name <- "oxy_2021"

library(indicspecies)
library(data.table)
# shouldn't need to load these libraries if using main.R
library(dplyr)
library(tidyr)
library(stringr)
library(glue)


# Create output dir if it doesn't exist yet
ind_sp_dir <- glue("results/{dataset_name}/indicator_species_analysis")
ifelse(!dir.exists(ind_sp_dir),
       dir.create(ind_sp_dir, mode = "777"), FALSE)

##################
## Regex expressions for different types
## 1) ABC{*}{e,u}_t{*}
##    - (\D\D\D(..)([e|u])_t(.)).*
## 2) ABC{*}_d{*}
##    - (\D\D\D(..)_d(.)).*
## 3) A_CTL{*}_d0 and A_BCD{*}_d{*}
##    - (\D_\D\D\D(..)_d(.)).*
## 4) ABC{*}_d{*}_t{*}
##    - (\D\D\D(..)_d(.)_t(.)).*

## split the string after first 3 characters.
## if there is an underscore in that then (3)
## else
##   count how many underscores are in it 
##   if there are two then (4) 
##   else
##     count how many letters are in it
##     if there are 5 letters then (1)
##     if there are 4 letters then (2)
##     else some error
##################

id_type <- function(tt){
  
  # 1) ABC{*}{e,u}_t{*}
  # 2) ABC{*}_d{*}
  # 3) A_CTL{*}_d0 and A_BCD{*}_d{*}
  # 4) ABC{*}_d{*}_t{*}
  
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
    if(how_many2 == 2){
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


regex_exp <- function(id_num, tt){
  #nt <- tt %>%
  #select(-c(name, taxRank, taxID, depth, taxLineage))
  d_nm <- toupper(substr(dataset_name, 1, nchar(dataset_name)-5))
  if(id_num == 1){
    ctl <- paste0(substr(d_nm, 1,1), "_CTL")
    reg <- paste0("((?:", d_nm, "|", ctl, ")(..)([e|u])_t(.)).*") # dont need to know whether it is D_flp or D_ctl because d0 will be if control
  }
  if(id_num == 2){
    reg <- paste0("(", d_nm, "(..)_d(.)).*")
  }
  if(id_num == 3){
    reg <- paste0("(", d_nm, "(..)_d(.)).*")
  }
  if(id_num == 4){
    reg <- paste0("(", d_nm, "(..)_d(.)_t(.)).*")
  }
  return(reg)
}


# write a function to identify type
# write a function to return regex expression
# write a function to return column names
# gonna need a way to dynamically do the mutate function

# calculates proportions and convert to data frame
calc_prop <- function(d) {
  sample_col <- select(d, c("sample", "replicate", "treatment"))
  prop_data <- as.matrix(select(d, -c("sample", "replicate", "treatment"))) %>%
    apply(MARGIN = 1, FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample, replicate, treatment)
  
  # calculate how proportion of sample attributed to each taxa
  
  return(prop_data)
}
# pass in filepath to tables file
run_indicator_analysis <- function(table_string, analysis_func, rank_symbol="all") {
  file <- paste0("results/",dataset_name,"/",table_string,"_counts.csv")
  table <- read.csv(file)
  # Pivot tables into wide format (1 column per taxa, 1 row per sample) and
  # create metadata category columns
  if(rank_symbol == "all"){
    tt <- table
  } else{
    tt <- filter(table, taxRank == rank_symbol)
  }
  d_nm <- substr(dataset_name, 1, nchar(dataset_name)-5)
  tt <- tt %>%
    pivot_longer(
      cols = contains(d_nm), # pivot columns containing "HBB"
      names_to = c("sample", "replicate", "treatment"), # make 3 new columns named these. These will contain parts of the names of tt cols
      names_pattern = "(OXY0(.)_d(.)).*", # split up tt column names to bin parts in names_to
      values_to = "read_count"
    ) %>%
    mutate(
      treatment = case_when(
        treatment == "0" ~ "Control",
        treatment == "1" ~ "Oxytetracycline"
      ),
      replicate = paste0("Rep ", replicate)#,
      #name = case_when(
        #name == "Actinobacteria" & taxRank == "C" ~ "Actinobacteria_class",
        #TRUE ~ name
      #) this part looks to be specific to the CTX datasets
    )
  # tt <- if(rank_symbol == "all") { tt }
  # else { filter(tt, taxRank == rank_symbol) }
  # browser()
  # Actinobacteria features twice in tt. once as taxrank C and once as taxrank P. We need to find a way to distinguish between those so we dont lose data
  name_and_id <- paste0(tt$name, " (", tt$taxID, ")")
  tt <- cbind(name_and_id, tt)
  
  tt <- tt %>% 
    select(-c(taxRank, taxLineage, taxID, depth, name))  %>%
    pivot_wider(names_from = "name_and_id", values_from = "read_count")
  
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
  
  # adjusted p-vals----------------
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
  
  # write to file
  sink(glue("{ind_sp_dir}/{func_out_str}_{table_string}_treatment_indicators_{rank_symbol}.txt"))
  summary(inv_treat)
  cat("\n\n\nAdjusted p-values\n")
  cat("-----------------\n")
  inv_treat.sign
  sink()
  sink(glue("{ind_sp_dir}/{func_out_str}_{table_string}_replicate_indicators_{rank_symbol}.txt"))
  summary(inv_rep)
  cat("\n\n\nAdjusted p-values\n")
  cat("-----------------\n")
  inv_rep.sign
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
