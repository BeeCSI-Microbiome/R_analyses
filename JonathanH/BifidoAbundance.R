library(dplyr)
library(tidyr)
library(ggplot2)

# read data
ctx_data <- read.delim(file = '2020_ctx_kraken2/ctx_kraken_genus_data.tsv',
                       header = TRUE,
                       sep = '\t')

# clean data
clean_ctx_data <- select(ctx_data, -taxRank, -taxID, -lineage) %>%
  filter(name %in% c("Gilliamella", 
                     "Frischella", 
                     "Snodgrassella", 
                     "Lactobacillus", 
                     "Bifidobacterium")) %>%
  pivot_longer(!name, names_to = "sample", values_to = "count") %>%
  pivot_wider(names_from = "name", values_from = "count")

# add treatment info
treatments <- c("Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI",
                "Control", "CLO", "THI")
clean_ctx_data$treatment <- treatments

# for ordering on plots
order <- c("Control", "CLO", "THI")

# adjust factor levels for ordering
clean_ctx_data$treatment <- factor(clean_ctx_data$treatment,
                                   levels = order)

# add replicate info
replicates <- c("Rep 2","Rep 2","Rep 2",
                "Rep 3","Rep 3","Rep 3",
                "Rep 4","Rep 4","Rep 4",
                "Rep 5","Rep 5","Rep 5",
                "Rep 6","Rep 6","Rep 6")
clean_ctx_data$replicate <- replicates

# reorder data frame
clean_ctx_data <- clean_ctx_data[, c(1,8,7,2,3,4,5,6)]

bifido_plot <- ggplot(clean_ctx_data, 
                      aes(x = replicate, 
                          y = Bifidobacterium, 
                          fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_discrete(limits = order) +
  labs(title = "Bifidobacterium Abundance",
       x = "Replicate",
       y = "Read Counts",
       fill = "Treatment")

bifido_plot

# save plot as svg
# svg("BifidoAbundance.svg")
# bifido_plot
# dev.off()