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
datapath <- 'results/thi_2020/plot_data/thi_raw_clade.csv'
treat_names <- c("Control","Acute", "Sublethal")
rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")
plot_title <- "Caged Control vs Acute THI - Genus DA"
plot_file_name <- 'thi_control_acute_genus_treatrep.png'
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

genus_obj <- prep_data(datapath, treat_names, rep_names, "G")

# run ancombc
mod <- ancombc(genus_obj,
               formula = 'treatment+replicate',
               p_adj_method = "BH")

# optionally inspect result using
# mod$res

# Visualize ---------------------------------------------------------------

# only get treatment related betas, adj p-vals, and diffs
df <- data.frame(mod$res) %>%
  select(contains('treatment')) %>%
  select(contains('beta')|contains('q_val')|contains('diff'))

# should work regardless of number of treatments and for every treatment
# paste('treatment', treat_names[2], sep = '')

df <- mutate(df, acute_diff_abun = ifelse(df$q_val.treatmentAcute < 0.05,
                                          ifelse(df$beta.treatmentAcute >= 0,
                                                 'Increase',
                                                 'Decrease'),
                                          'No Change'))

df <- mutate(df, sublethal_diff_abun = ifelse(df$q_val.treatmentSublethal < 0.05,
                                              ifelse(df$beta.treatmentSublethal >= 0,
                                                     'Increase',
                                                     'Decrease'),
                                              'No Change'))

# adjust factor levels for plotting
df$acute_diff_abun <- factor(df$acute_diff_abun,
                             levels = c('Increase', 'No Change', 'Decrease'))
df$sublethal_diff_abun <- factor(df$sublethal_diff_abun,
                                 levels = c('Increase', 'No Change', 'Decrease'))

da_plot <- ggplot(df,
                  aes(x = beta.treatmentAcute,
                      y = -log10(q_val.treatmentAcute),
                      colour = acute_diff_abun)) +
  geom_point(alpha=0.5, size=3) +
  # may need to modify colour values based on whether there are any increases/decreases
  scale_color_manual(values=c('blue', 'black', 'red')) +
  # should check that limits allow all points to be visible
  xlim(c(-2,2)) +
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05), lty=4,col="black",lwd=0.8) +
  labs(x = "ln(fold change)",
       y = "-log10(adj. p-value)",
       title = plot_title) +
  theme(legend.position = "right",
        legend.title = element_blank())

da_plot


# Save Results ------------------------------------------------------------
ggsave(da_plot,
       filename = plot_file_name)

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

