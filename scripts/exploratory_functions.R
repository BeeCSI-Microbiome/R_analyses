# Function scripts
# Author(s): Jonathan Ho, Lance Lansing

# Load Packages -----------------------------------------------------------
import('dplyr')
import('tidyr')
import('ggplot2')
import('vegan')
import('stats', 'aggregate')
import("glue")

# Global Variables --------------------------------------------------------
# Taxa of interest - common
interest_list <- c('Lactobacillus Firm-4',
                   'Lactobacillus Firm-5',
                   'Other Lactobacillus',
                   'Gilliamella apicola',
                   'Gilliamella apis',
                   'Pantoea agglomerans',
                   'Bifidobacterium',
                   'Snodgrassella alvi',
                   'Frischella perrara',
                   'Bartonella apis',
                   'Melissococcus plutonius',
                   'Paenibacillus larvae',
                   'Bacteria')

# Add taxa specific to current dataset
interest_list <- append(interest_list, c('Gilliamella apis'))

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

# Reorders table using taxa column
order_taxa <- function(d) {
  taxa_order <- get_taxa_order(d)
  d$taxa <- factor(d$taxa, levels=taxa_order)
  d
}

# Returns a list of taxa of interest in order of descending average percent
# but with 'Other Bacteria' always last
get_taxa_order <- function(d) {
  pdlong <- filter(d, taxa!='Other Bacteria')
  
  pd_avg <- aggregate(pdlong[,c('value')], list(pdlong$taxa), mean) %>%
    arrange(value)
  
  taxa_order <- append(pd_avg$Group.1, 'Other Bacteria')
}


# Relative Abundance ------------------------------------------------------

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

# Returns relative abundance for taxa of interest
export("make_interest_abundance")
make_interest_abundance <- function(data, treat_names, rep_names, dataset_name) {
  plot_data <- filter(data, name %in% interest_list) %>%
    add_other_bac() %>%
    tidy_data() %>%
    calc_prop() %>%
    treat_reps(treat_names, rep_names)
    
  plot <- tidy_to_long(plot_data) %>%
    order_taxa() %>%
    plot_interest_abundance()
  
  ggsave(plot = plot,
         filename = glue('results/{dataset_name}/interest_abundance.png'),
         bg = 'white')
  utils::write.csv(plot_data,
                   file = glue('results/{dataset_name}/plot_data/interest_taxa_proportions.csv'),
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
    order_taxa() %>%
    plot_interest_abundance()
  
  thi_plot <- filter(process, treatment %in% thi_treat) %>%
    tidy_to_long() %>%
    order_taxa() %>%
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
make_all_alpha_plots <- function(data, treat_names, rep_names, dataset_name) {
  plot_data <- prep_alpha_data(data, treat_names, rep_names)
  
  utils::write.csv(plot_data,
                   file = glue('results/{dataset_name}/plot_data/alpha_div.csv'),
                   row.names = F)
  
  plot_1 <- plot_alpha(plot_data, "Shannon")
  plot_2 <- plot_alpha(plot_data, "Simpson")
  plot_3 <- plot_alpha(plot_data, "Inv_Simpson")
  plot_4 <- plot_alpha(plot_data, "Evenness")
  
  ggsave(plot = plot_1, filename = glue('results/{dataset_name}/alpha_div_shannon.png'), bg = 'white')
  ggsave(plot = plot_2, filename = glue('results/{dataset_name}/alpha_div_simpson.png'), bg = 'white')
  ggsave(plot = plot_3, filename = glue('results/{dataset_name}/alpha_div_inv_simp.png'), bg = 'white')
  ggsave(plot = plot_4, filename = glue('results/{dataset_name}/alpha_div_evenness.png'), bg = 'white')
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

# plots the nmds for activity 1 data
plot_nmds_1 <- function(data, h_var, plot_title) {
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

# plots the nmds for activity 2 data
plot_nmds_2 <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  
  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data = hull, alpha = 0.5) +
    geom_point(aes(shape = replicate), size = 3) +
    aes(fill = !!hull_var) +
    labs(title = plot_title,
         x = "NMDS1", 
         y = "NMDS2",
         shape = "Replicate",
         fill = tools::toTitleCase(h_var))
  
  plot
}

# plots and saves nmds' with separate hulls for replicate and treatment
act_1_nmds <- function(genus_data, speci_data, dataset_name) {
  genus_treat_nmds <- plot_nmds_1(genus_data,
                                "treatment",
                                "Genus NMDS - Brays Curtis")
  genus_reps_nmds <- plot_nmds_1(genus_data,
                               "replicate",
                               "Genus NMDS - Brays Curtis")
  speci_treat_nmds <- plot_nmds_1(speci_data,
                                "treatment",
                                "Species NMDS - Bray Curtis")
  speci_reps_nmds <- plot_nmds_1(speci_data,
                               "replicate",
                               "Species NMDS - Bray Curtis")
  
  ggsave(plot = genus_treat_nmds, 
         filename = glue('results/{dataset-name}/genus_treat_nmds.png'), bg = 'white')
  ggsave(plot = genus_reps_nmds, 
         filename = glue('results/{dataset-name}/genus_reps_nmds.png'), bg = 'white')
  ggsave(plot = speci_treat_nmds, 
         filename = glue('results/{dataset-name}/speci_treat_nmds.png'), bg = 'white')
  ggsave(plot = speci_reps_nmds, 
         filename = glue('results/{dataset-name}/speci_reps_nmds.png'), bg = 'white')
}

# plots and saves nmds' for treatment only
act_2_nmds <- function(genus_data, speci_data, dataset_name) {
  genus_treat_nmds <- plot_nmds_2(genus_data,
                                "treatment",
                                "Genus NMDS - Brays Curtis")
  speci_treat_nmds <- plot_nmds_2(speci_data,
                                "treatment",
                                "Species NMDS - Bray Curtis")
  ggsave(plot = genus_treat_nmds,
         filename = glue('results/{dataset_name}/genus_treat_nmds.png'), bg = 'white')
  ggsave(plot = speci_treat_nmds,
         filename = glue('results/{dataset_name}/speci_treat_nmds.png'), bg = 'white')
}

# make and save nmds plots for genus and species levels using read data
# uses raw reads to calculate proportions
export("make_nmds_plots")
make_nmds_plots <- function(data, treat_names, rep_names, dataset_name) {
  genus_data <- filter(data, taxRank == "G") %>%
    prep_nmds_data(treat_names, rep_names)
  speci_data <- filter(data, taxRank == "S") %>%
    prep_nmds_data(treat_names, rep_names)
  
  utils::write.csv(genus_data,
                   file = glue('results/{dataset_name}/plot_data/genus_nmds.csv'),
                   row.names = F)
  utils::write.csv(speci_data,
                   file = glue('results/{dataset_name}/plot_data/species_nmds.csv'),
                   row.names = F)
  
  # split act 1 and 2 nmds here, based on number of treatments
  if (length(treat_names) > 2) {
    # can make convex hulls for both replicate and treatment
    act_1_nmds(genus_data, speci_data, dataset_name)
  } else {
    # can only make hulls for treatment
    act_2_nmds(genus_data, speci_data, dataset_name)
  }
  
}

# testing nmds using only interest list
export("interest_nmds")
interest_nmds <- function(data, treat_names, rep_names, dataset_name) {
  interest_data <- filter(data, name %in% interest_list) %>%
    prep_nmds_data(treat_names, rep_names)
  
  utils::write.csv(interest_data,
                   file = glue('results/{dataset_name}/plot_data/interest_nmds.csv'),
                   row.names = F)
  
  int_treat_nmds <- plot_nmds_1(interest_data,
                              "treatment",
                              "Interest NMDS - Bray Curtis")
  int_rep_nmds <- plot_nmds_1(interest_data,
                            "replicate",
                            "Interest NMDS - Bray Curtis")
  
  ggsave(plot = int_treat_nmds, 
         filename = glue('results/{dataset_name}/interest_treat_nmds.png'), bg = 'white')
  ggsave(plot = int_rep_nmds, 
         filename = glue('results/{dataset_name}/interest_rep_nmds.png'), bg = 'white')
}

# testing nmds using all taxa data
export("all_taxa_nmds")
all_taxa_nmds <- function(data, treat_names, rep_names, dataset_name) {
  all_data <- data %>%
    prep_nmds_data(treat_names, rep_names)
  
  utils::write.csv(all_data,
                   file = glue('results/{dataset_name}/plot_data/all_taxa_nmds.csv'),
                   row.names = F)
  
  all_treat_nmds <- plot_nmds_1(all_data,
                              "treatment",
                              "All NMDS - Bray Curtis")
  all_rep_nmds <- plot_nmds_1(all_data,
                            "replicate",
                            "All NMDS - Bray Curtis")
  
  ggsave(plot = all_treat_nmds, filename = glue('results/{dataset_name}/all_treat_nmds.png'), bg = 'white')
  ggsave(plot = all_rep_nmds, filename = glue('results/{dataset_name}/all_rep_nmds.png'), bg = 'white')
}