library(dplyr)
library(tidyr)
library(ggplot2)

# read data
ctx_data <- read.delim(file = 'ctx_kraken_genus_precent_data.tsv',
                       header = TRUE,
                       sep = '\t')

# clean data
bifido_data <- select(ctx_data, -taxRank, -taxID, -Max, -lineage) %>%
  filter(name == "Bifidobacterium") %>%
  pivot_longer(!name, names_to = "sample", values_to = "count") %>%
  pivot_wider(names_from = "name", values_from = "count")

# add treatment info
treatments <- c("Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI")
bifido_data$treatment <- treatments

# for ordering on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
bifido_data$treatment <- factor(bifido_data$treatment,
                                   levels = order)

# add replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")
bifido_data$replicate <- replicates

bifido_plot <- ggplot(bifido_data, 
                      aes(x = replicate, 
                          y = Bifidobacterium, 
                          fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(limits = order) +
  labs(title = "Bifidobacterium Percent Abundance",
       x = "Replicate",
       y = "Read Counts",
       fill = "Treatment")

bifido_plot

# save plot as svg
# svg("Bifido_Percent_Abundance.svg")
# bifido_plot
# dev.off()
