# Function scripts
# Author(s): Jonathan Ho, Lance Lansing

# Load Packages -----------------------------------------------------------
import('dplyr')
import('tidyr')
import('ggplot2')
import('vegan')

# Global Variables --------------------------------------------------------
interest_list <- c('Lactobacillus Firm-4',
                   'Lactobacillus Firm-5',
                   'Other Lactobacillus',
                   'Gilliamella apicola',
                   'Bifidobacterium',
                   'Snodgrassella alvi',
                   'Frischella perrara',
                   'Bartonella apis',
                   'Melissococcus plutonius',
                   'Paenibacillus larvae',
                   'Bacteria')

# Wrangling Functions -----------------------------------------------------
# cleans data into tidy format
tidy_data <- function(d) {
  clean_data <- select(d, -taxRank, -lineage) %>%
    pivot_longer(!name, names_to = "sample", values_to = "value") %>%
    pivot_wider(names_from = "name", values_from = "value")
    
  return(clean_data)
}

# calculates and adds Other Bacteria col
add_other_bac <- function(d) {
  column_index <- grepl('Reads', names(d))
  # Get the sum counts for taxa of interest, except Bacteria (domain), then 
  # subtract their sum from Bacteria clade count. This leaves the Bacteria row
  # representing all taxa of non-interest 
  df <- d[d$name!='Bacteria', column_index]
  d[d$name=='Bacteria', column_index] <-
    d[d$name=='Bacteria', column_index] - as.list(colSums(df, na.rm=T))
  # Change the grouping name
  d[d$name=='Bacteria',]$name <- 'Other Bacteria'
  
  return(d)
}


# calculates proportions and convert to data frame
calc_prop <- function(d) {
  sample_col <- select(d, "sample")
  prop_data <- select(d, -"sample") %>%
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
    t() %>%
    as.data.frame() %>%
    mutate(sample_col)
  
  return(prop_data)
}

# add others column and select for core genera
filter_core <- function(d) {
  core_data <- mutate(d, 
                      Others = 100 - (Gilliamella +
                                        Snodgrassella + 
                                        Bifidobacterium +
                                        Lactobacillus +
                                        Frischella)) %>%
    select("Gilliamella",
           "Snodgrassella",
           "Bifidobacterium",
           "Lactobacillus",
           "Frischella",
           "Others",
           "sample")
  return(core_data)
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

# convert tidy data into long for abundance plot setup
tidy_to_long <- function(d) {
  long_data <- pivot_longer(d,
                            cols = !c("sample", "treatment", "replicate"),
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
  
  plot_data <- filter(data, taxRank == "G") %>%
    tidy_data() %>%
    calc_prop() %>%
    filter_core() %>%
    treat_reps(treat_names, rep_names)
    
  
  plot <- tidy_to_long(plot_data) %>%
    plot_genera_abundance(plot_data)
  
  ggsave(plot = plot, filename = 'results/relative_abundance.png', bg = 'white')
  utils::write.csv(plot_data,
                   file = 'results/plot_data/genera_proportions.csv',
                   row.names = F)
}

# Returns relative abundance for taxa of interest
export("make_interest_abundance")
make_interest_abundance <- function(data, treat_names, rep_names) {
  plot_data <- filter(data, name %in% interest_list) %>%
    add_other_bac() %>%
    tidy_data() %>%
    calc_prop() %>%
    treat_reps(treat_names, rep_names)
  
  plot <- tidy_to_long(plot_data) %>%
    plot_interest_abundance()
  
  ggsave(plot = plot, filename = 'results/interest_abundance.png', bg = 'white')
  utils::write.csv(plot_data,
                   file = 'results/plot_data/interest_taxa_proportions.csv',
                   row.names = F)
}

# separate bar plots for each treatment vs control
# CTX experiment specific, not necessary for any other experiment
export("make_separate_ctx_bars")
make_separate_ctx_bars <- function(data, treat_names, rep_names) {
  clo_treat <- c("Control", "CLO")
  thi_treat <- c("Control", "THI")
  
  process <- filter(data, name %in% interest_list) %>%
    add_other_bac() %>%
    tidy_data() %>%
    calc_prop() %>%
    treat_reps(treat_names, rep_names) 
  
  clo_plot <- filter(process, treatment %in% clo_treat) %>%
    tidy_to_long() %>%
    plot_interest_abundance()
  
  thi_plot <- filter(process, treatment %in% thi_treat) %>%
    tidy_to_long() %>%
    plot_interest_abundance()
  
  ggsave(plot = clo_plot, filename = 'results/clo_abundance.png', bg = 'white')
  ggsave(plot = thi_plot, filename = 'results/thi_abundance.png', bg = 'white')
}

# Alpha Diversity ---------------------------------------------------------
# calc alpha metrics
# assumes data in tidy format and no extra columns
calc_diversity_df <- function(d){
  sample_col <- select(d, "sample")
  sample_data <- select(d, -"sample")
  observed_richness <- specnumber(sample_data)
  invsimpson <- diversity(sample_data, index="invsimpson")
  simpson <- diversity(sample_data, index="simpson")
  shannon <- diversity(sample_data, index="shannon")
  evenness <- shannon/log(observed_richness)
  div_df <- data.frame(
    ID = sample_col,
    Observed_Richness = observed_richness,
    Inv_Simpson = invsimpson,
    Simpson = simpson,
    Shannon = shannon,
    Evenness = evenness
  )
  return(div_df)
}

# preps and saves alpha div data
prep_alpha_data <- function(d, treat_names, rep_names) {
  plot_data <- filter(d, taxRank == "S") %>%
    tidy_data() %>%
    calc_diversity_df() %>%
    treat_reps(treat_names, rep_names)
  
  utils::write.csv(plot_data,
                   file = 'results/plot_data/alpha_div.csv',
                   row.names = F)
  
  return(plot_data)
}

# plot alpha diversity
plot_alpha <- function(d, alpha) {
  index <- sym(alpha)
  alpha_plot <- ggplot(d, aes(x = treatment, y = !!index)) +
    geom_boxplot() +
    labs(title = "Alpha Diversity",
         x = "Treatment")
}

# Makes all alpha div bar plots and saves them in results
export("make_all_alpha_plots")
make_all_alpha_plots <- function(data, treat_names, rep_names) {
  plot_data <- prep_alpha_data(data, treat_names, rep_names)
  
  plot_1 <- plot_alpha(plot_data, "Shannon")
  plot_2 <- plot_alpha(plot_data, "Simpson")
  plot_3 <- plot_alpha(plot_data, "Inv_Simpson")
  plot_4 <- plot_alpha(plot_data, "Evenness")
  
  ggsave(plot = plot_1, filename = 'results/alpha_div_shannon.png', bg = 'white')
  ggsave(plot = plot_2, filename = 'results/alpha_div_simpson.png', bg = 'white')
  ggsave(plot = plot_3, filename = 'results/alpha_div_inv_simp.png', bg = 'white')
  ggsave(plot = plot_4, filename = 'results/alpha_div_evenness.png', bg = 'white')
}


# Beta Diversity ----------------------------------------------------------
# calculates nmds using Bray Curtis dissimilarity using only taxa data
calc_nmds <- function(data) {
  set.seed(1)
  sample_col <- select(data, "sample")
  
  nmds_data <- select(data, -"sample") %>%
    as.matrix() %>%
    metaMDS(distance = "bray") %>%
    scores() %>%
    as.data.frame()%>%
    mutate(sample_col)
  
  return(nmds_data)
}

# preps nmds data
prep_nmds_data <- function(d, treat_names, rep_names) {
  genus_data <- tidy_data(d) %>%
    calc_prop() %>%
    calc_nmds() %>%
    treat_reps(treat_names, rep_names)

  return(genus_data)
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

# make and save nmds plots for genus and species levels using read data
# uses raw reads to calculate proportions
export("make_nmds_plots")
make_nmds_plots <- function(data, treat_names, rep_names) {
  genus_data <- filter(data, taxRank == "G") %>%
    prep_nmds_data(treat_names, rep_names)
  speci_data <- filter(data, taxRank == "S") %>%
    prep_nmds_data(treat_names, rep_names)
  
  utils::write.csv(genus_data,
                   file = 'results/plot_data/genus_nmds.csv',
                   row.names = F)
  utils::write.csv(speci_data,
                   file = 'results/plot_data/species_nmds.csv',
                   row.names = F)
  
  genus_treat_nmds <- plot_nmds(genus_data,
                                "treatment",
                                "Genus NMDS - Brays Curtis")
  genus_reps_nmds <- plot_nmds(genus_data,
                               "replicate",
                               "Genus NMDS - Brays Curtis")
  speci_treat_nmds <- plot_nmds(speci_data,
                                "treatment",
                                "Species NMDS - Bray Curtis")
  speci_reps_nmds <- plot_nmds(speci_data,
                               "replicate",
                               "Species NMDS - Bray Curtis")
  
  ggsave(plot = genus_treat_nmds, filename = 'results/genus_treat_nmds.png', bg = 'white')
  ggsave(plot = genus_reps_nmds, filename = 'results/genus_reps_nmds.png', bg = 'white')
  ggsave(plot = speci_treat_nmds, filename = 'results/speci_treat_nmds.png', bg = 'white')
  ggsave(plot = speci_reps_nmds, filename = 'results/speci_reps_nmds.png', bg = 'white')
}
