# Function scripts
# Author(s): Jonathan Ho

# Load Packages -----------------------------------------------------------
packs <- c('dplyr', 'tidyr', 'ggplot2', 'vegan')
lapply(packs, library, character.only=T)

# Wrangling Functions -----------------------------------------------------
# cleans data into tidy format for genus or species level
tidy_data <- function(d, r) {
  if(r=="G") {
    to_remove <- c("Apis")
  } else if(r=="S") {
    to_remove <- c("Apis mellifera")
  }
  clean_data <- filter(d, taxRank == r) %>%
    select(-taxRank, -taxID, -Max, -lineage) %>%
    pivot_longer(!name, names_to = "sample", values_to = "value") %>%
    filter(!name %in% to_remove) %>%
    pivot_wider(names_from = "name", values_from = "value")%>%
  return(clean_data)
}

# calculates proportions and convert to data frame
# ignores the first col (typically sample names)
calc_prop <- function(d){
  prop_data <- apply(d[,2:ncol(d)],
                       MARGIN = 1,
                       FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame()
  return(prop_data)
}

# select for core taxa and add Others col
filter_core <- function(d) {
  core_data <- select(d,
                      "Gilliamella",
                      "Snodgrassella",
                      "Bifidobacterium",
                      "Lactobacillus",
                      "Frischella") %>%
    mutate(Others = 100 - (Gilliamella +
                           Snodgrassella + 
                           Bifidobacterium +
                           Lactobacillus +
                           Frischella))
}

# add in treatment and replicate cols
# assumes same number of replicates for each treatment and vice versa. 
# assumes that samples are grouped by replicates first, then treatments.
# Example below:
# Rep 1 TreatmentA, Rep 1 TreatmentB, Rep 1 TreatmentC, 
# Rep 2 TreatmentA, Rep 2 TreatmentB, Rep 2 TreatmentC, ...
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

# convert tidy data into long for abundance plot setup
# assumes last 2 cols are treatment and replicate, so ignored
tidy_to_long <- function(d) {
  long_data <- pivot_longer(d,
                            cols = 1:(ncol(d)-2),
                            names_to = "taxa",
                            values_to = "value")
  return(long_data)
}



# Relative Abundance ------------------------------------------------------
# plots relative abundance data at genus level
# uses the taxa, value, treatment, and replicate column from data
plot_genera_abundance <- function(d) {
  abundance_plot <- ggplot(d, 
                           aes(x = treatment,
                               y = value,
                               fill = factor(taxa,
                                             levels = c("Bifidobacterium",
                                                        "Frischella",
                                                        "Gilliamella",
                                                        "Lactobacillus",
                                                        "Snodgrassella",
                                                        "Others")))) +
    geom_bar(stat = "identity", colour = "black") +
    facet_grid(~replicate) +
    labs(title = "Relative Abundance",
         x = "Treatment",
         y = "Percentage (%)",
         fill = "Genus") +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1))
}

# Returns a relative abundance plot at the genus level
# see specific functions for details and assumptions
make_genera_abundance <- function(datapath, treat_names, rep_names) {
  data <- read.delim(file = datapath,
                     header = TRUE,
                     sep = '\t')
  plot <- tidy_data(data, "G") %>%
    calc_prop() %>%
    filter_core() %>%
    treat_reps(treat_names, rep_names) %>%
    tidy_to_long() %>%
    plot_genera_abundance()
    
  plot
}

# Relative Abundance Tests ------------------------------------------------
# necessary info
# datapath <- '../2020_ctx_kraken2/ctx_kraken_all_rawread_uncollapsed.tsv'
# treat_names <- c("Control", "CLO", "THI")
# rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")

# data <- read.delim(file = datapath,
#                    header = TRUE,
#                    sep = '\t')

# tests for genera abundance
# test_1 <- tidy_data(data, "G")
# test_2 <- calc_prop(test_1)
# test_3 <- filter_core(test_2)
# test_4 <- treat_reps(test_3, treat_names, rep_names)
# test_5 <- tidy_to_long(test_4)
# test_6 <- plot_genera_abundance(test_5, plot_title)
# test_6

# make_genera_abundance(datapath,
#                       treat_names,
#                       rep_names)


# Alpha Diversity ---------------------------------------------------------
# calc alpha metrics
# assumes data in tidy format and no extra columns
calc_diversity_df <- function(x){
  observed_richness <- specnumber(x[,2:ncol(x)])
  invsimpson <- diversity(x[,2:ncol(x)], index="invsimpson")
  simpson <- diversity(x[,2:ncol(x)], index="simpson")
  shannon <- diversity(x[,2:ncol(x)], index="shannon")
  evenness <- shannon/log(observed_richness)
  div_df <- data.frame(
    ID = x$sample,
    Observed_Richness = observed_richness,
    Inv_Simpson = invsimpson,
    Simpson = simpson,
    Shannon = shannon,
    Evenness = evenness
  )
  div_df
}

# plot alpha diversity
plot_alpha <- function(d, alpha) {
  index <- sym(alpha)
  alpha_plot <- ggplot(d, aes(x = treatment, y = !!index)) +
    geom_boxplot() +
    labs(title = "Alpha Diversity Plot",
         x = "Treatment")
  
  alpha_plot
}

# Returns an alpha diversity plot using specific index
# see specific functions for details and assumptions
make_alpha_bars <- function(datapath, treat_names, rep_names, alpha_index) {
  data <- read.delim(file = datapath,
                     header = TRUE,
                     sep = '\t')
  plot <- tidy_data(data, "S") %>%
    calc_diversity_df() %>%
    treat_reps(treat_names, rep_names) %>%
    plot_alpha(alpha_index)
  
  plot
}

# Alpha Diversity Tests ---------------------------------------------------
# necessary info
# datapath <- '../2020_ctx_kraken2/ctx_kraken_all_rawread_uncollapsed.tsv'
# treat_names <- c("Control", "CLO", "THI")
# rep_names <- c("Rep 2", "Rep 3", "Rep 4", "Rep 5", "Rep 6")
# plot_title <- "CTX Abundance Using Percent(%) Data"

# data <- read.delim(file = datapath,
#                    header = TRUE,
#                    sep = '\t')

# tests for alpha diversity
# test_1 <- tidy_data(data, "S")
# test_2 <- calc_diversity_df(test_1)
# test_3 <- treat_reps(test_2, treat_names, rep_names)
# test_4 <- plot_alpha(test_3, "Shannon")
# test_4

# make_alpha_bars(datapath,
#                 treat_names,
#                 rep_names,
#                 "Inv_Simpson")
