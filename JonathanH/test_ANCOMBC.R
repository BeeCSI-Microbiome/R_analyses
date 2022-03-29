# This script is for testing ANCOM_BC as another method for differential
# abundance. This is to be done before adding to the main R workflow.
# Takes some functions from exploratory_functions.R

# Author(s): Jonathan Ho
# Updated: Mar 29, 2022

library(tidyverse)
library(ANCOMBC)
library(microbiome)

# User Defined Variables --------------------------------------------------
datapath <- 'results/clo_2020/plot_data/clo_raw_clade.csv'
treat_names <- c("Control","Acute", "Sublethal")
rep_names <- c("Rep 1", "Rep 2", "Rep 3", "Rep 5", "Rep 6")
lab_names <- c("Lab 1", "Lab 2")
plot_title <- "Caged Control vs Acute CLO - Species DA"
plot_file_name <- 'clo_control_acute_species_treatrep.png'
summary_file_name <- "clo_2020_species_mod2_trunc.csv"

# Functions ---------------------------------------------------------------
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

# adds lab metadata for ctx experiment
add_labs <- function(d, lab_names, treat_names) {
  # calculate num of each lab from number of treatments given
  num_lab1 <- length(treat_names)*3
  num_lab2 <- length(treat_names)*2
  
  # add lab1
  lab = c(rep(lab_names[1], num_lab1))
  
  # add lab2
  lab = c(lab, rep(lab_names[2], num_lab2))

  d$lab <- lab
  return(d)
}


# Setup -------------------------------------------------------------------

data <- read_csv(datapath) %>%
  filter(taxRank == "S")

abun_data <- select(data, -taxRank, -lineage) %>%
  column_to_rownames('name') %>%
  as.matrix()

metadata <- tidy_data(data) %>%
  select(sample) %>%
  treat_reps(treat_names, rep_names) %>%
  column_to_rownames('sample') %>%
  add_labs(lab_names, treat_names)

OTU <- otu_table(abun_data, taxa_are_rows = T)
samples <- sample_data(metadata)

phylo_obj <- phyloseq(OTU, samples)


# Call --------------------------------------------------------------------
# 3 different models differing in their formulas
# mod1 <- ancombc(phylo_obj,
#                 formula = 'treatment+replicate+lab',
#                 p_adj_method = "BH",
#                 zero_cut = 0.9,
#                 lib_cut = 0,
#                 group = 'treatment',
#                 struc_zero = F,
#                 neg_lb = F,
#                 tol = 1e-05,
#                 max_iter = 100,
#                 conserve = F,
#                 alpha = 0.05,
#                 global = T)


# includes replicate as confounding variable, uses rep 2 as baseline
mod2 <- ancombc(phylo_obj,
                formula = 'treatment+replicate',
                p_adj_method = "BH",
                group = 'treatment',
                global = T)

# includes lab as confounding variable, uses lab 1 as baseline
# mod3 <- ancombc(phylo_obj,
#                 formula = 'treatment+lab',
#                 p_adj_method = "BH",
#                 group = 'treatment',
#                 global = T)

# # save full output
# write.csv(mod2$res,
#           file = "ctx_2020_clo_species_mod2.csv",
#           row.names = T)
# 
# write.csv(mod3$res,
#           file = "ctx_2020_clo_species_mod3.csv",
#           row.names = T)


# mod1$res
# mod1$res_global

# mod2$res$diff_abn
# mod2$res_global
# 
# mod3$res
# mod3$res_global


# Visualize ---------------------------------------------------------------

df <- data.frame(mod2$res) %>%
  select(contains('treatmentAcute')|contains('treatmentSublethal')) %>%
  select(contains('beta')|contains('q_val')|contains('diff'))

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
#              ln_fold_change= beta.treatmentTHI) %>%
#   select(-diff_abn.treatmentTHI)
# 
# df <- df[with(df, order(adj_p_val, -ln_fold_change)),]

write.csv(df,
          file = summary_file_name,
          row.names = T)
