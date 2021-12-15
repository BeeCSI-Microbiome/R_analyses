# Function scripts
# Author(s): Jonathan Ho

# Load Packages -----------------------------------------------------------
import('dplyr')
import('tidyr')
import('ggplot2')
import('vegan')

# Wrangling Functions -----------------------------------------------------
# cleans data into tidy format
tidy_data <- function(d) {
  clean_data <- select(d, -taxRank, -lineage) %>%
    pivot_longer(!name, names_to = "sample", values_to = "value") %>%
    pivot_wider(names_from = "name", values_from = "value")
    
  clean_data
}

# calculates proportions and convert to data frame
# subsets out the first col (typically sample names)
calc_prop <- function(d){
  prop_data <- apply(d[,2:ncol(d)],
                     MARGIN = 1,
                     FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame()
  
  prop_data
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

# plots relative abundance data for taxa of interest
# uses the taxa, value, treatment, and replicate column from data
plot_interest_abundance <- function(d) {
  abundance_plot <- ggplot(d, 
                           aes(x = treatment,
                               y = value,
                               fill = taxa)) +
    geom_bar(stat = "identity", colour = "black") +
    facet_grid(~replicate) +
    labs(title = "Relative Abundance",
         x = "Treatment",
         y = "Percentage (%)",
         fill = "Taxa") +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1))
}

# Returns a relative abundance plot at the genus level
# see specific functions for details and assumptions
export('make_genera_abundance')
make_genera_abundance <- function(data, treat_names, rep_names) {
  plot <- filter(data, taxRank == "G") %>%
    tidy_data() %>%
    calc_prop() %>%
    filter_core() %>%
    treat_reps(treat_names, rep_names) %>%
    tidy_to_long() %>%
    plot_genera_abundance()
  
  plot
  ggsave(plot = plot, filename = 'results/relative_abundance.png', bg = 'white')
}

# Returns relative abundance for taxa of interest
export("make_interest_abundance")
make_interest_abundance <- function(data, treat_names, rep_names, interest_list) {
  plot <- filter(data, name %in% interest_list) %>%
    tidy_data() %>%
    calc_prop() %>%
    treat_reps(treat_names, rep_names) %>%
    tidy_to_long() %>%
    plot_interest_abundance()
  
  plot
  ggsave(plot = plot, filename = 'results/interest_abundance.png', bg = 'white')
  
}

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
    labs(title = "Alpha Diversity",
         x = "Treatment")
  
  alpha_plot
}

# Returns an alpha diversity plot using specific index
# see specific functions for details and assumptions
make_alpha_bars <- function(data, treat_names, rep_names, alpha_index) {
  plot <- filter(data, taxRank == "S") %>%
    tidy_data() %>%
    calc_diversity_df() %>%
    treat_reps(treat_names, rep_names) %>%
    plot_alpha(alpha_index)
  
  plot
}

# Makes all alpha div bar plots and saves them in results
export("make_all_alpha_plots")
make_all_alpha_plots <- function(data, treat_names, rep_names) {
  plot_1 <- make_alpha_bars(data, treat_names, rep_names, "Shannon")
  plot_2 <- make_alpha_bars(data, treat_names, rep_names, "Simpson")
  plot_3 <- make_alpha_bars(data, treat_names, rep_names, "Inv_Simpson")
  plot_4 <- make_alpha_bars(data, treat_names, rep_names, "Evenness")
  ggsave(plot = plot_1, filename = 'results/alpha_div_shannon.png', bg = 'white')
  ggsave(plot = plot_2, filename = 'results/alpha_div_simpson.png', bg = 'white')
  ggsave(plot = plot_3, filename = 'results/alpha_div_inv_simp.png', bg = 'white')
  ggsave(plot = plot_4, filename = 'results/alpha_div_evenness.png', bg = 'white')
}


# Beta Diversity ----------------------------------------------------------
# calculates nmds using Bray Curtis distances using only taxa data
calc_nmds <- function(data) {
  set.seed(1)
  nmds_data <- data %>%
    as.matrix() %>%
    metaMDS(distance = "bray") %>%
    scores() %>%
    as.data.frame()
  
  nmds_data
}

# plots the nmds
plot_nmds <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  
  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(shape = 21, size = 3) +
    geom_polygon(data = hull, alpha = 0.5) +
    aes(fill = !!hull_var) +
    labs(title = plot_title,
         x = "NMDS1", 
         y = "NMDS2",
         fill = tools::toTitleCase(h_var))
  
  plot
}

# makes an nmds plot using tidy raw read data
make_nmds <- function(data, treat_names, rep_names, h_var, plot_title) {
  
  # calc_prop subsets data to remove samples column
  nmds_data <- calc_prop(data) %>%
    calc_nmds()
  
  # add back samples column
  nmds_data$sample <- data$sample
  
  plot <- nmds_data %>%
    treat_reps(treat_names, rep_names) %>%
    plot_nmds(h_var, plot_title)
  
  plot
}


# make and save nmds plots for genus and species levels using read data
# uses raw reads to calculate proportions (for now)
export("make_nmds_plots")
make_nmds_plots <- function(data, treat_names, rep_names) {
  genus_data <- filter(data, taxRank == "G") %>%
    tidy_data()
  speci_data <- filter(data, taxRank == "S") %>%
    tidy_data()
  
  genus_treat_nmds <- make_nmds(genus_data,
                                treat_names,
                                rep_names,
                                "treatment",
                                "Genus NMDS - Brays Curtis")
  genus_reps_nmds <- make_nmds(genus_data,
                               treat_names,
                               rep_names,
                               "replicate",
                               "Genus NMDS - Brays Curtis")
  speci_treat_nmds <- make_nmds(speci_data,
                                treat_names,
                                rep_names,
                                "treatment",
                                "Species NMDS - Bray Curtis")
  speci_reps_nmds <- make_nmds(speci_data,
                               treat_names,
                               rep_names,
                               "replicate",
                               "Species NMDS - Bray Curtis")
  
  ggsave(plot = genus_treat_nmds, filename = 'results/genus_treat_nmds.png', bg = 'white')
  ggsave(plot = genus_reps_nmds, filename = 'results/genus_reps_nmds.png', bg = 'white')
  ggsave(plot = speci_treat_nmds, filename = 'results/speci_treat_nmds.png', bg = 'white')
  ggsave(plot = speci_reps_nmds, filename = 'results/speci_reps_nmds.png', bg = 'white')
}
