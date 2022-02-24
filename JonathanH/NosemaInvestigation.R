# Script for investigating Nosema in data samples
# TODO:'s show recommended fields that should be changed for each analysis

library(dplyr)
library(tidyr)
library(ggplot2)
library(glue)

# TODO: change file path
dataset_name <- 'ctx_2020'
datapath <- 'results/ctx_2020/ctx_2020_read_summary.csv'
treat_names <- c('Control', 'CLO', 'THI')
rep_names <- c('Rep 1', 'Rep 2', 'Rep 3', 'Rep 4', 'Rep 5')

# function for adding treatment and replicate info
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

# read data
data <- read.delim(file = datapath,
                       header = TRUE,
                       sep = ',')

# clean data
data <- select(data, sample, contains('Nosema')) %>%
  mutate(percent_classified_Nosema = percent_classified_Nosema_apis + percent_classified_Nosema_ceranae)
data <- treat_reps(data, treat_names, rep_names)

# plot data
read_plot <- ggplot(data,
                     aes(x = replicate,
                         y = percent_classified_Nosema,
                         fill = treatment))+
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Percent Classified Nosema Reads",
       x = "Replicate",
       y = "Reads",
       fill = "Treatment")

read_plot

ggsave(plot = read_plot, filename = glue('results/{dataset_name}/nosema.png'), bg = 'white')

# wilcox test (nonparametric since we don't know distribution)
control_vs_clo <- filter(data, treatment %in% c("Control", "CLO"))
wilcox.test(percent_classified_Nosema~treatment, data=control_vs_clo)

