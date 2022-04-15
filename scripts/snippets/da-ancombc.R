# Script for running ANCOM-BC for for pairwise differential 
# abundance analysis and producing volcano plots.
# Takes two functions from exploratory_functions.R for data wrangling:
# tidy_data and treat_reps
# Takes as input: raw clade counts (should not be scaled/normalized)
# Uses treatment+replicate for as the model, but this can be 
# customized in the ancombc function itself.

# Author(s): Jonathan Ho
# Updated: Apr 11, 2022

library(tidyverse)
library(ANCOMBC)
library(microbiome)

# User Defined Variables --------------------------------------------------
datapath <- 'results/clo_2020/plot_data/clo_raw_clade.csv'
treat_names <- c("Control","Acute","Sublethal")
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 5", "Rep 6")
base_title <- 'Caged Control vs'
dataset_name <- 'clo_2020'
summary_file_name <- "thi_2020_genus_mod_trunc.csv"


# Functions from exploratory_functions.R ----------------------------------

# cleans data into tidy format
tidy_data <- function(d) {
  clean_data <- select(d, -taxRank, -lineage) %>%
    pivot_longer(!name, names_to = "sample", values_to = "value") %>%
    pivot_wider(names_from = "name", values_from = "value")
  
  return(clean_data)
}

# add in treatment and replicate cols
# assumes same number of replicates for each treatment and vice versa. 
# assumes that samples are grouped by replicates first, then treatments.
# Example below:
# Rep 1 TreatmentA, Rep 1 TreatmentB, Rep 1 TreatmentC, 
# Rep 2 TreatmentA, Rep 2 TreatmentB, Rep 2 TreatmentC, ...
# Can double check that these are correct by comparing with samples col
treat_reps <- function(d, treat_names, rep_names) {
  num_treats <- length(treat_names)
  num_reps <- length(rep_names)
  
  treatments <- rep(treat_names, num_reps)
  replicates <- c()
  for (r in 1:num_reps) {
    replicates = c(replicates, rep(rep_names[r], num_treats))
  }
  
  d$treatment <- treatments
  d$replicate <- replicates
  
  # adjust factor levels for future plotting
  d$treatment <- factor(d$treatment,
                        levels = treat_names)
  
  return(d)
}


# ANCOM-BC ----------------------------------------------------------------
prep_data <- function(datapath, treat_names, rep_names, rank) {
  data <- read_csv(datapath) %>%
    filter(taxRank == rank)
  
  # setup abundance data
  abun_data <- select(data, -taxRank, -lineage) %>%
    column_to_rownames('name') %>%
    as.matrix()
  
  # setup metadata
  metadata <- tidy_data(data) %>%
    select(sample) %>%
    treat_reps(treat_names, rep_names) %>%
    column_to_rownames('sample')
  
  # combine into phyloseq object
  OTU <- otu_table(abun_data, taxa_are_rows = T)
  samples <- sample_data(metadata)
  
  phylo_obj <- phyloseq(OTU, samples)
  
  return(phylo_obj)
  
}



# optionally inspect result using
# mod$res

# Visualize ---------------------------------------------------------------
# setup data frame for plotting and later saving
visualize_save <- function(ancombc_mod, treat_names) {
  
  # only get treatment related betas, adj p-vals, and diffs
  df <- data.frame(mod$res) %>%
    select(contains('treatment')) %>%
    select(contains('beta')|contains('q_val')|contains('diff'))
  
  # loop through treatment names and do everything inside
  controls <- c('Control', 'unexposed')
  just_treats <- treat_names[!treat_names %in% controls]
  for (i in just_treats) {
    # setup identifiable strings for i-th treatment
    diff_abun_i <- paste('diff_abun_', i, sep = '') %>%
      sym()
    q_val_i <- paste('q_val.treatment', i, sep = '') %>%
      sym()
    beta_i <- paste('beta.treatment', i, sep = '') %>%
      sym()
    plot_title_i <- paste(base_title, i, '-', dataset_name)
    
    # determines whether DA are increases or decreases
    df <- mutate(df, !!diff_abun_i := ifelse(df[,grep(q_val_i, 
                                                      colnames(df))] < 0.05,
                                             ifelse(df[,grep(beta_i,
                                                             colnames(df))] >= 0,
                                                    'Increase',
                                                    'Decrease'),
                                             'No Change'))
    
    # adjust factor levels for plotting
    df[,grep(diff_abun_i, colnames(df))] <- factor(df[,grep(diff_abun_i, colnames(df))],
                                                   levels = c('Increase', 'No Change', 'Decrease'))
    
    # plot and save
    da_plot <- volc_plot(df, beta_i, q_val_i, diff_abun_i, plot_title_i)
    
    ggsave(da_plot,
           filename = paste(plot_title_i, '.png', sep = ''))
  }

}
          
# plots volcano plot for given data and column names in data frame
volc_plot <- function(df, beta_i, q_val_i, diff_abun_i, plot_title_i) {

  # determine x range
  x_default = 1.25
  x_data = max(abs(df[,grep(beta_i,colnames(df))]))
  x_use = max(x_default, x_data)
  
  # determine y range
  y_default = 1.5
  y_data = -log10(min(df[,grep(q_val_i,colnames(df))]))
  y_use = max(y_default, y_data) 
  
  # plot
  da_plot <- ggplot(df,
                    aes(x = !!beta_i,
                        y = -log10(!!q_val_i),
                        colour = !!diff_abun_i)) +
    geom_point(alpha=0.5, size=3) +
    scale_color_manual(values=c('Increase' = 'blue',
                                'No Change' = 'black',
                                'Decrease' = 'red')) +
    xlim(0-x_use, 
         0+x_use) +
    ylim(0, y_use) +
    geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
    geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +
    labs(x = "ln(fold change)",
         y = "-log10(adj. p-value)",
         title = plot_title_i) +
    theme(legend.position = "right",
          legend.title = element_blank())
  
  return(da_plot)
}

# TODO: use single function to combine these 3 and do both genus and species
genus_obj <- prep_data(datapath, treat_names, rep_names, "G")

# run ancombc
mod <- ancombc(genus_obj,
               formula = 'treatment+replicate',
               p_adj_method = "BH")
visualize_save(mod, treat_names)



######## Testing

# only get treatment related betas, adj p-vals, and diffs
df <- data.frame(mod$res) %>%
  select(contains('treatment')) %>%
  select(contains('beta')|contains('q_val')|contains('diff'))

# loop through treatment names and do everything inside
controls <- c('Control', 'unexposed')
just_treats <- treat_names[!treat_names %in% controls]

# setup identifiable strings for i-th treatment
diff_abun_i <- paste('diff_abun_', just_treats[2], sep = '') %>%
  sym()
q_val_i <- paste('q_val.treatment', just_treats[2], sep = '') %>%
  sym()
beta_i <- paste('beta.treatment', just_treats[2], sep = '') %>%
  sym()
plot_title_i <- paste(base_title, just_treats[2], '-', dataset_name)

# determines whether DA are increases or decreases
df <- mutate(df, !!diff_abun_i := ifelse(df[,grep(q_val_i, 
                                                  colnames(df))] < 0.05,
                                         ifelse(df[,grep(beta_i,
                                                         colnames(df))] >= 0,
                                                'Increase',
                                                'Decrease'),
                                         'No Change'))

# adjust factor levels for plotting
df[,grep(diff_abun_i, colnames(df))] <- factor(df[,grep(diff_abun_i, colnames(df))],
                                               levels = c('Increase', 'No Change', 'Decrease'))

# determine x range
x_default = 1.25
x_data = max(abs(df[,grep(beta_i,colnames(df))]))
x_use = max(x_default, x_data)

# determine y range
y_default = 1.5
y_data = -log10(min(df[,grep(q_val_i,colnames(df))]))
y_use = max(y_default, y_data) 

ggplot(df,
       aes(x = !!beta_i,
           y = -log10(!!q_val_i),
           colour = !!diff_abun_i)) +
  geom_point(alpha=0.5, size=3) +
  scale_color_manual(values=c('Increase' = 'blue',
                              'No Change' = 'black',
                              'Decrease' = 'red')) +
  xlim(0-x_use, 
       0+x_use) +
  ylim(0, y_use) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +
  labs(x = "ln(fold change)",
       y = "-log10(adj. p-value)",
       title = plot_title_i) +
  theme(legend.position = "right",
        legend.title = element_blank())

# sort and reorder data frame
df <- rename(df,
             adj_p_val = q_val.treatmentTHI ,
             ln_fold_change = beta.treatmentTHI) %>%
  select(-diff_abn.treatmentTHI)

df <- df[with(df, order(q_val_i, -ln_fold_change)),]

# Save Results ------------------------------------------------------------
# TODO: incorporate saving df and sorting by significant q values

# sort and reorder data frame
# df <- rename(df,
#              adj_p_val = q_val.treatmentTHI ,
#              ln_fold_change = beta.treatmentTHI) %>%
#   select(-diff_abn.treatmentTHI)
# 
# df <- df[with(df, order(adj_p_val, -ln_fold_change)),]

write.csv(df,
          file = summary_file_name,
          row.names = T)

