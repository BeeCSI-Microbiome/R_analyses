# Function scripts
# Author(s): Jonathan Ho, Lance Lansing

# Load Packages -----------------------------------------------------------
import("dplyr")
import("tidyr")
import("ggplot2")
import("vegan")
import("stats", "aggregate")
import("glue")
import("stringr")

# Scope Variables --------------------------------------------------------
# Taxa of interest - common
#'*is there an equivalent list for ARGs?  I think no*
# interest_list <- c("Lactobacillus Firm-4",
#                    "Lactobacillus Firm-5",
#                    "Other Lactobacillus",
#                    "Gilliamella apicola",
#                    "Gilliamella apis",
#                    "Bifidobacterium",
#                    "Snodgrassella alvi",
#                    "Frischella perrara",
#                    "Bartonella apis",
#                    "Melissococcus plutonius",
#                    "Paenibacillus larvae",
#                    "Bacteria")


# Wrangling Functions -----------------------------------------------------
# cleans data into tidy format
#'*commented out to preserve original code*'
#'#'*select() is removing those columns, which we don't have in the first place*
# tidy_data <- function(d) {
#   clean_data <- select(d, -taxRank, -taxLineage, -taxID, -depth) %>%
#     pivot_longer(!name, names_to = "sample", values_to = "value") %>%
#     pivot_wider(names_from = "name", values_from = "value")
# 
#   return(clean_data)
# }
#'*updated. see screenshot 463 to see result*
tidy_data_samples <- function(data) {
  clean_data <- select(data, -geneID, -type, -mechanism, -gene) %>%
    # browser()
    pivot_longer(!class, names_to = "sample", values_to = "value") 
  # browser()
  clean_data[is.na(clean_data)] <- 0
  clean_data = group_by(clean_data, sample, class) %>% 
    summarise(value=sum(value)) %>%
    # browser()
    pivot_wider(names_from = "class", values_from = "value")  
  browser()

    # data.frame(crop = c("HBB"))
  clean_data = ungroup(clean_data)
    return(clean_data)
}

# tidy_data_samples <- function(data) {
#   clean_data <- select(data, -geneID, -type, -mechanism, -gene) %>%
#     # browser()
#     pivot_longer(!class, names_to = "sample", values_to = "value") %>% 
#     # browser()
#     group_by(class, sample) %>% 
#     summarise(value=sum(value)) %>%
#     # browser()
#     pivot_wider(names_from = "class", values_from = "value") 
#   browser()
#   
#   # data.frame(crop = c("HBB"))
#   
#   return(clean_data)
# }

tidy_data_all <- function(data) {
  clean_data <- select(data, -geneID, -type, -mechanism, -gene) %>%
    # browser()
    pivot_longer(!class, names_to = "sample", values_to = "value") 
    # browser()
  clean_data[is.na(clean_data)] <- 0
  clean_data = group_by(clean_data, sample, class) %>% 
    summarise(value=sum(value)) %>%
    # browser()
    pivot_wider(names_from = "class", values_from = "value") 
  #browser()
  crop_label = substr(clean_data$sample,1,3)
  year_label = str_sub(clean_data$sample, -2, -1)
  #browser()
  clean_data$crop = paste(crop_label,year_label, sep = "")
  
  #browser()
  #'*does this distort the values?*
  # clean_data[is.na(clean_data)] <- 0
  # browser()
  # data.frame(crop = c("HBB"))
  clean_data = ungroup(clean_data)
  return(clean_data)
}

tidy_data_province <- function(data) {
  clean_data <- select(data, -geneID, -type, -mechanism, -gene) %>%
    # browser()
    pivot_longer(!class, names_to = "sample", values_to = "value") 
    # browser()
  clean_data[is.na(clean_data)] <- 0
  clean_data = group_by(clean_data, sample, class) %>% 
    summarise(value=sum(value)) %>%
    # browser()
    pivot_wider(names_from = "class", values_from = "value") 
  browser()
  crop_label = substr(clean_data$sample,1,3)
  year_label = str_sub(clean_data$sample, -2, -1)
  treatment_label = substr(clean_data$sample,6,6)
  # browser()
  clean_data$crop = paste(crop_label,year_label,treatment_label, sep = "")
  clean_data$province = str_sub(clean_data$sample, -5, -4)
  # browser()
  #'*does this distort the values?*
  # clean_data[is.na(clean_data)] <- 0
  # browser()
  # data.frame(crop = c("HBB"))
  clean_data = ungroup(clean_data)
  return(clean_data)
}

tidy_data_crop_av <- function(data) {
  clean_data <- select(data, -geneID, -type, -mechanism, -gene) %>%  
    pivot_longer(!class, names_to = "sample", values_to = "value") %>%    
    group_by(class) %>% 
    summarise(value=sum(value)) %>% 
    pivot_wider(names_from = "class", values_from = "value") %>% 
    data.frame(crop = c("HBB"))
  
  return(clean_data)
}

#'*don't need this - not dealing with taxa*
# calculates and adds Other Bacteria col
# add_other_bac <- function(d) {
#   column_index <- !(names(d) %in% c("taxRank","taxLineage","taxID","depth","name"))
#   # Get the sum counts for taxa of interest, except Bacteria (domain), then 
#   # subtract their sum from Bacteria clade count. This leaves the Bacteria row
#   # representing all taxa of non-interest 
#   df <- d[d$name!="Bacteria", column_index]
#   d[d$name=="Bacteria", column_index] <-
#     d[d$name=="Bacteria", column_index] - as.list(colSums(df, na.rm=T))
#   # Change the grouping name
#   d[d$name=="Bacteria",]$name <- "Other Bacteria"
#   
#   return(d)
# }

# calculates and adds Other Bacteria col
add_other_bac <- function(d) {
  column_index <- !(names(d) %in% c("taxRank","taxLineage","taxID","depth","name"))
  # Get the sum counts for taxa of interest, except Bacteria (domain), then
  # subtract their sum from Bacteria clade count. This leaves the Bacteria row
  # representing all taxa of non-interest
  df <- d[d$name!="Bacteria", column_index]
  d[d$name=="Bacteria", column_index] <-
    d[d$name=="Bacteria", column_index] - as.list(colSums(df, na.rm=T))
  # Change the grouping name
  d[d$name=="Bacteria",]$name <- "Other ARGs"

  return(d)
}



# calculates proportions and convert to data frame
# calc_prop <- function(d) {
#   sample_col <- select(d, "sample")
#   prop_data <- select(d, -"sample") %>%
#     apply(MARGIN = 1,
#           FUN = function(x) x / sum(x) * 100) %>%
#     t() %>%
#     as.data.frame() %>%
#     mutate(sample_col) %>%
#     relocate(sample)
#   
#   return(prop_data)
# }

calc_prop_samples <- function(d) {
  sample_col <- select(d, "sample")
  plot_data_samples <- select(d, -"sample") %>% 
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
  # browser()
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample) 

  #utils::write.csv(plot_data_samples, file = "car-speeds-cleaned.csv")
  # browser()
  return(plot_data_samples)
}

calc_prop_all <- function(d) {
  # crop_col <- select(d, "crop")
  sample_col <- select(d, "sample")
  crop_col <- select(d, "crop")
  #browser()
  plot_data_samples <- select(d, -"sample", -"crop") %>% 
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
    #browser()
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample) 
  mutate(crop_col) %>%
    relocate(crop)
  #browser()
  #utils::write.csv(plot_data_samples, file = "car-speeds-cleaned.csv")
  # browser()
  return(plot_data_samples)
}

calc_prop_province <- function(d) {
  # crop_col <- select(d, "crop")
  sample_col <- select(d, "sample")
  crop_col <- select(d, "crop")
  province_col = select(d, "province")
  browser()
  plot_data_samples <- select(d, -"sample", -"crop", -"province") %>% 
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
    #browser()
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample) 
    mutate(crop_col) %>%
    relocate(crop)
    mutate(province_col) %>%
      relocate(province)
  #browser()
  #utils::write.csv(plot_data_samples, file = "car-speeds-cleaned.csv")
  # browser()
  return(plot_data_samples)
}

calc_prop_crop_av <- function(d) {
  crop_col <- select(d, "crop")
  plot_data_crop_av <- select(d, -"crop") %>%
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>% 
    t() %>%
    as.data.frame() %>% 
    mutate(crop_col) %>%
    relocate(crop)
  
  return(plot_data_crop_av)
}
# add in treatment and replicate cols
# Samples gathered from column names are treated with regular expressions to 
# extract replicate number and treatments.
province_treat_reps <- function(d, treatment_key) {
  # browser()
  #'*can probably get rid of treat_reps for all crop analysis*
  d <- d %>% mutate(replicate = extract_replicate_string(sample),
                    treatment = extract_treatment_string(sample, treatment_key))
  #'*is the following going to cause problems for the individual crop analyses?*
  crop_label = substr(d$sample,1,3)
  year_label = str_sub(d$sample, -2, -1)
  treatment_label = substr(d$sample,6,6)
  d$crop = paste(crop_label,year_label,treatment_label, sep = "")
  d$province = str_sub(d$sample, -5, -4)
  
  #browser()
  return(d)
}

all_treat_reps <- function(d, treatment_key) {
  # browser()
  #'*can probably get rid of treat_reps for all crop analysis*
  d <- d %>% mutate(replicate = extract_replicate_string(sample),
                    treatment = extract_treatment_string(sample, treatment_key))
  #'*is the following going to cause problems for the individual crop analyses?*
  crop_label = substr(d$sample,1,3)
  year_label = str_sub(d$sample, -2, -1)
  d$crop = paste(crop_label,year_label, sep = "")
  
  #browser()
  return(d)
}

treat_reps <- function(d, treatment_key) {
  # browser()
  #'*can probably get rid of treat_reps for all crop analysis*
  d <- d %>% mutate(replicate = extract_replicate_string(sample),
                    treatment = extract_treatment_string(sample, treatment_key))

  
  #browser()
  return(d)
}

extract_replicate_string <- function(sample_string) {
  digits <- str_extract(sample_string, "(?<=\\w{3,5})\\d\\d") %>% 
    str_remove("^0+")
    paste0("Rep ", digits)
}


extract_treatment_string <- function(sample_string, treatment_key) {
  if (all(str_detect(sample_string, "_d[[:alnum:]]+"))) {
    # Does sample string match activity 1 pattern?
    treatment_codes <-
      str_extract(sample_string, "(?<=_)d[[:alnum:]]+")
  } else if (all(str_detect(sample_string, "(?<=[[:upper:]]{3}\\d\\d)(e|u)"))) {
    # Or activity 2 pattern?
    treatment_codes <-
      str_extract(sample_string, "(?<=[[:upper:]]{3}\\d\\d)(e|u)")
  } else {
    stop("The sample column names do not match the implemented regex patterns")
  }
  
  if (all(treatment_codes %in% names(treatment_key))) {
    unlist(treatment_key[treatment_codes], use.names = FALSE)
  } else {
    missing_codes <-
      unique(treatment_codes)[!unique(treatment_codes) %in% names(treatment_key)]
    stop(paste(
          "The following treatment code(s) detected from samples names were not found in treatment key you provinceided:",
          paste(missing_codes, collapse = ", ")))
  }
}

# convert tidy data into long for abundance plot setup
# tidy_to_long <- function(d) {
#   long_data <- pivot_longer(d,
#                             cols = !c("sample", "treatment", "replicate"),
#                             names_to = "taxa",
#                             values_to = "value")
#   return(long_data)
# }


tidy_to_long_samples <- function(d) {
  long_data_samples <- pivot_longer(d,
                            cols = !c("sample", "treatment", "replicate", "crop"), #"crop"
                            names_to = "ARGs",
                            values_to = "value")
  browser()
  return(long_data_samples)
}

tidy_to_long_all <- function(d) {
  long_data_samples <- pivot_longer(d,
                                    cols = !c("sample", "treatment", "replicate", "crop"), #"crop"
                                    names_to = "ARGs",
                                    values_to = "value")
  return(long_data_samples)
}

tidy_to_long_crop_av <- function(d) {
  long_data_crop_av <- pivot_longer(d,
                                    cols = !c("crop"),
                                    names_to = "ARGs",
                                    values_to = "value")
  return(long_data_crop_av)
}

# Reorders table using taxa column
order_taxa <- function(d) {
  browser()
  taxa_order <- get_taxa_order(d)
  d$crop <- factor(d$crop, levels=taxa_order)
  
}

# Returns a list of taxa of interest in order of descending average percent
# but with "Other Bacteria" always last
# get_taxa_order <- function(d) {
#   pdlong <- filter(d, taxa!="Other Bacteria")
#   
#   pd_avg <- aggregate(pdlong[,c("value")], list(pdlong$taxa), mean) %>%
#     arrange(value)
#   
#   taxa_order <- append(pd_avg$Group.1, "Other Bacteria")
# }


get_taxa_order <- function(d) {
  #pdlong <- filter(d, taxa!="Other Bacteria")
  browser()
  d$crop = substr(d$sample,1,3)
  browser()
  pd_avg <- aggregate(d[,c("value")], list(d$crop), mean) %>%
    arrange(value)
  browser()
  #taxa_order <- append(pd_avg$Group.1, "Other Bacteria")
}


# Relative Abundance ------------------------------------------------------

# plots relative abundance data for taxa of interest
# uses the taxa, value, treatment, and replicate column from data
# plot_interest_abundance <- function(d) {
#   abundance_plot <- ggplot(d, 
#                            aes(x = treatment,
#                                y = value,
#                                fill = taxa)) +
#     geom_bar(stat = "identity", colour = "black") +
#     facet_grid(~replicate) +
#     labs(title = "Relative Abundance",
#          x = "Treatment",
#          y = "Percentage (%)",
#          fill = "Taxa") +
#     theme(axis.text.x = element_text(angle = 45,
#                                      vjust = 1,
#                                      hjust = 1))
# }

#'*creates the graphs with all the samples laid out. Works, except the graph isn't pretty.  Bars too thin, etc.*
plot_interest_abundance_samples <- function(d) {
  abundance_plot <- ggplot(d, 
                           aes(x = treatment,
                               y = value,
                               fill = ARGs)) +
    geom_bar(width = 1, stat = "identity", colour = "black") +
    facet_grid(~replicate) +
    labs(title = "Relative Abundance",
         x = "Treatment",
         y = "Percentage (%)",
         fill = "ARG Class") +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1))
}

#'*creates the graphs with all the samples averaged into one bar.  Needs work*
plot_interest_abundance_crop_av <- function(d) {
  abundance_plot <- ggplot(d, 
                           aes(x = crop,
                               y = value,
                               fill = ARGs)) +
    geom_bar(width = 1, stat = "identity", colour = "black") +
    labs(title = "Relative Abundance",
         x = "Crop",
         y = "Percentage (%)",
         fill = "ARG Class") +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1))
}

# Returns relative abundance for taxa of interest
#'*starting point, that calls all other functions in this script*
#'*"data" is being fed "tables[["raw_clade"]]" from main.r*
export("make_interest_abundance")
make_interest_abundance <- function(data, treatment_key, dataset_name, outdir) {
  # if (any(!is.na(additional_taxa))) {
  #   interest_list <- append(interest_list, additional_taxa)
  # }
  
  #'*original plot_data before my edits:*
  # plot_data <- filter(data, name %in% interest_list) %>%
  #   add_other_bac() %>%
  #   tidy_data() %>%
  #   calc_prop() %>%
  #   treat_reps(treatment_key)
  
  plot_data_samples <- tidy_data_samples(data) %>% 
    calc_prop_samples() %>% 
    treat_reps(treatment_key)
 
  # plot_data_crop_av <- tidy_data_crop_av(data) %>%
  #   calc_prop_crop_av()
  #treat_reps(treatment_key)

  #' #'*sends the relative abundance data up to the ggplot function*
  plot_samples <- tidy_to_long_samples(plot_data_samples) 
  browser()
  plot_samples <- order_taxa(plot_samples)
  browser()
  plot_samples <- plot_interest_abundance_samples(plot_samples)
  
  # plot_crop_av = tidy_to_long_crop_av(plot_data_crop_av) %>%
  #   order_taxa() %>% 
  #   plot_interest_abundance_crop_av()
  

  ggsave(plot = plot_samples,
         filename = glue("{outdir}/relative_abundance_samples_plot.png"),
         bg = "white")
  
  # ggsave(plot = plot_crop_av,
  #        filename = glue("{outdir}/relative_abundance_crop_av_plot.png"),
  #        bg = "white")
  utils::write.csv(plot_data_samples,
                   file = glue("{outdir}/relative_abundance_proportions_samples.csv"),
                   row.names = F)
  plot_samples
  # utils::write.csv(plot_data_crop_av,
  #                  file = glue("{outdir}/relative_abundance_proportions_crop_av.csv"),
  #                  row.names = F)
}

# separate bar plots for each treatment vs control
# CTX experiment specific, not necessary for any other experiment
export("make_separate_ctx_bars")
make_separate_ctx_bars <- function(data, treatment_key, dataset_name,
                                   additional_taxa, outdir) {
  
  interest_list <- append(interest_list, additional_taxa)
  
  clo_treat <- c("Control", "CLO")
  thi_treat <- c("Control", "THI")
  
  process <- filter(data, name %in% interest_list) %>%
    add_other_bac() %>%
    tidy_data() %>%
    calc_prop() %>%
    treat_reps(treatment_key) 
  
  clo_plot <- filter(process, treatment %in% clo_treat) %>%
    tidy_to_long() %>%
    order_taxa() %>%
    plot_interest_abundance()
  
  thi_plot <- filter(process, treatment %in% thi_treat) %>%
    tidy_to_long() %>%
    order_taxa() %>%
    plot_interest_abundance()
  
  ggsave(plot = clo_plot, filename = glue("{outdir}/clo_abundance.png"), bg = "white")
  ggsave(plot = thi_plot, filename = glue("{outdir}/thi_abundance.png"), bg = "white")
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

# preps alpha div data
prep_alpha_data <- function(d, treatment_key) {
  plot_data <- filter(d, taxRank == "S") %>%
    tidy_data() %>%
    calc_diversity_df() %>%
    treat_reps(treatment_key)
}

# runs and saves kruskal wallis on shannon index for both 
# replicates and treatments
alpha_KW_test <- function(d, dataset_name, index, outdir) {
  heading <- glue("{index} Index - Kruskal Wallis Results:\n")
  rep_alpha <- stats::kruskal.test(stats::formula(glue("{index}~replicate")), d)
  treat_alpha <- stats::kruskal.test(stats::formula(glue("{index}~treatment")), d)
  
  cat(heading, file = glue("{outdir}/alpha_stats.txt"), append = T)
  utils::capture.output(rep_alpha, file = glue("{outdir}/alpha_stats.txt"), append = T)
  utils::capture.output(treat_alpha, file = glue("{outdir}/alpha_stats.txt"), append = T)
}

# plot alpha diversity
plot_alpha <- function(d, alpha, dataset_name) {
  index <- sym(alpha)
  alpha_plot <- ggplot(d, aes(x = treatment, y = !!index)) +
    geom_boxplot() +
    labs(title = "Alpha Diversity",
         x = "Treatment")
}

# Makes all alpha div bar plots, calculates simple stats and saves them in results
export("make_all_alpha_plots")
make_all_alpha_plots <- function(data, treatment_key, dataset_name, outdir) {
  plot_data <- prep_alpha_data(data, treatment_key)
  
  utils::write.csv(plot_data,
                   file = glue("{outdir}/alpha_div_data.csv"),
                   row.names = F)
  
  file.create(glue("{outdir}/alpha_stats.txt"))
  alpha_KW_test(plot_data, dataset_name, "Shannon", outdir)
  alpha_KW_test(plot_data, dataset_name, "Simpson", outdir)

  plot_1 <- plot_alpha(plot_data, "Shannon", dataset_name)
  plot_2 <- plot_alpha(plot_data, "Simpson", dataset_name)
  plot_3 <- plot_alpha(plot_data, "Inv_Simpson", dataset_name)
  plot_4 <- plot_alpha(plot_data, "Evenness", dataset_name)
  
  ggsave(plot = plot_1, filename = glue("{outdir}/alpha_div_shannon.png"), bg = "white")
  ggsave(plot = plot_2, filename = glue("{outdir}/alpha_div_simpson.png"), bg = "white")
  ggsave(plot = plot_3, filename = glue("{outdir}/alpha_div_inverse_simpson.png"), bg = "white")
  ggsave(plot = plot_4, filename = glue("{outdir}/alpha_div_evenness.png"), bg = "white")
}


# Beta Diversity ----------------------------------------------------------
# calculates nmds using Bray Curtis dissimilarity using only taxa data
calc_nmds <- function(data) {
  set.seed(1)
  sample_col <- select(data, "sample")
  nmds_data <- select(data, -"sample") %>% #, -"crop"
    as.matrix() %>%
    metaMDS(distance = "bray") %>%
    scores()
  nmds_data <- nmds_data$sites %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample)
  # browser()
  return(nmds_data)
}

#'*not sure if this is doing what it is supposed to - should it still be removing sample?*
#'* running into a bug with the NAs in the dataframe, I think*
all_calc_nmds <- function(data) {
  #'*Setting a seed in R means to initialize a pseudorandom number generator*
  browser()
  set.seed(1)
  sample_col <- select(data, "sample")
  nmds_data <- select(data, -"sample") %>%  #, -"crop"?
  # browser()
  #'*as.matrix converts a data.table into a matrix, optionally using one of the columns in the data.table as the matrix row names*
    as.matrix() %>%
    metaMDS(distance = "bray") %>%
  #'*Function to access either species or site scores for specified axes in some ordination methods*
    scores()
  #browser()
  nmds_data <- nmds_data$sites %>% 
  #browser()  
  as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample)
  #browser()
  return(nmds_data)
}

pair_calc_nmds <- function(data) {
  #'*Setting a seed in R means to initialize a pseudorandom number generator*
  browser()
  set.seed(1)
  sample_col <- select(data, "sample")
  crop_col = select(data, "crop")
  prov_col = select(data, "crop")
    nmds_data <- select(data, -"sample", -"crop", -"province") %>%  #, -"crop"?
    # browser()
    #'*as.matrix converts a data.table into a matrix, optionally using one of the columns in the data.table as the matrix row names*
    as.matrix() %>%
    metaMDS(distance = "bray") %>%
    #'*Function to access either species or site scores for specified axes in some ordination methods*
    scores()
  #browser()
  nmds_data <- nmds_data$sites %>% 
    #browser()  
    as.data.frame() %>%
    mutate(sample_col) %>%
    mutate(crop_col) %>%
    mutate(prov_col) %>%
    relocate(sample)
  #browser()
  return(nmds_data)
}

# runs anosim and saves results in a text file in results folder
#'*taxa_level is just the name of the graph, ie, "Crop's Samples"*
calc_ano <- function(d, group_data, taxa_level, dataset_name, outdir) {
  heading <- paste(taxa_level, "ANOSIM results:\n")
  ano_data <- select(d, -"sample", -"crop") %>%
    as.matrix()
  rep_ano <- anosim(ano_data,
                    group_data$replicate,
                    distance = "bray",
                    permutations = 9999)
  
  treat_ano <- anosim(ano_data,
                      group_data$treatment,
                      distance = "bray",
                      permutations = 9999)
  
  cat(heading,
      file = glue("{outdir}/anosim.txt"), append = T)
  utils::capture.output(
    rep_ano,
    file = glue("{outdir}/anosim.txt"), append = T)
  utils::capture.output(
    treat_ano,
    file = glue("{outdir}/anosim.txt"), append = T)
}

#'*group_data$crop contains both crop name and year, so the anosim IS comparing correctly*
all_calc_ano <- function(d, group_data, taxa_level, dataset_name, outdir) {
  heading <- paste(taxa_level, "ANOSIM results:\n")
  ano_data <- select(d, -"sample", -"crop") %>%
    as.matrix()
  crop_ano <- anosim(ano_data,
                    group_data$crop,
                    distance = "bray",
                    permutations = 9999)
  # browser()
  cat(heading,
      file = glue("{outdir}/anosim.txt"), append = T)
  utils::capture.output(
    crop_ano,
    file = glue("{outdir}/anosim.txt"), append = T)
}

province_calc_ano <- function(d, group_data, taxa_level, dataset_name, outdir) {
  heading <- paste(taxa_level, "ANOSIM results:\n")
  #browser()
  ano_data <- select(d, -"sample", -"crop", -"province") %>%
    as.matrix()
  province_ano <- anosim(ano_data,
                     group_data$province,
                     distance = "bray",
                     permutations = 9999)
  # browser()
  cat(heading,
      file = glue("{outdir}/anosim.txt"), append = T)
  utils::capture.output(
    province_ano,
    file = glue("{outdir}/anosim.txt"), append = T)
}

# preps nmds data and runs ANOSIM on data
#'*taxa_level is just the name of the graph, ie, "Crop's Samples"*
prep_and_ano <- function(d, treatment_key, taxa_level, dataset_name, outdir) {
  prop_data_samples <- tidy_data_samples(d) %>%
    calc_prop_samples()
  
  # browser()
  
  group_data <- prop_data_samples %>%
    calc_nmds() %>%
    treat_reps(treatment_key)
  # browser()
 
  prop_data_samples$crop = substr(prop_data_samples$sample,1,3)
  # browser()
  calc_ano(prop_data_samples, group_data, taxa_level, dataset_name, outdir)
  
  return(group_data)
}

all_prep_and_ano <- function(d, treatment_key, all_crops_title, dataset_name, outdir) {
  all_prop_data_samples <- tidy_data_all(d) %>%
    calc_prop_all()

  # browser()
  
  group_data <- all_prop_data_samples %>%
    all_calc_nmds() %>%
    all_treat_reps(treatment_key)
  # browser()
  
  all_prop_data_samples$crop = substr(all_prop_data_samples$sample,1,3)
  # browser()
  all_calc_ano(all_prop_data_samples, group_data, all_crops_title, dataset_name, outdir)
  
  return(group_data)
}

province_prep_and_ano <- function(d, treatment_key, province_title, dataset_name, outdir) {
  browser()
  #'*at this point HBB21 sample names have the format HBB05u_t2>, while all other crops are CORe_t2_M_ON_20.  ie, have lost the year and province info*
  province_prop_data_samples <- tidy_data_province(d) 
  browser()
  province_prop_data_samples = calc_prop_province(province_prop_data_samples)
  
  browser()
  
  group_data <- province_prop_data_samples %>%
    all_calc_nmds() %>%
    province_treat_reps(treatment_key)
  # browser()
  
  province_prop_data_samples$crop = substr(province_prop_data_samples$sample,1,3)
  province_prop_data_samples$province = substr(province_prop_data_samples$sample,-5,-4)
  # browser()
  province_calc_ano(province_prop_data_samples, group_data, province_title, dataset_name, outdir)
  
  return(group_data)
}


BC_AB_pair_prep_and_ano <- function(d, treatment_key, province_title, dataset_name, outdir) {
  browser()
  
  pair_province_prop_data_samples <- tidy_data_province(d) 
  browser()
  pair_province_prop_data_samples = calc_prop_province(pair_province_prop_data_samples)
  
  browser()
  
  pair_province_prop_data_samples$crop = substr(pair_province_prop_data_samples$sample,1,3)
  pair_province_prop_data_samples$province = str_sub(pair_province_prop_data_samples$sample,-5,-4)
  browser()
  BC_AB = pair_province_prop_data_samples %>% 
    filter(str_detect(province, "BC|AB", negate = FALSE))
  browser()
  
  group_data <- BC_AB %>%
    pair_calc_nmds() %>%
    province_treat_reps(treatment_key)
  # browser()
  province_calc_ano(BC_AB, group_data, province_title, dataset_name, outdir)
  
  return(group_data)
}

QC_AB_pair_prep_and_ano <- function(d, treatment_key, province_title, dataset_name, outdir) {
  browser()
  
  pair_province_prop_data_samples <- tidy_data_province(d) 
  browser()
  pair_province_prop_data_samples = calc_prop_province(pair_province_prop_data_samples)
  
  browser()
  
  pair_province_prop_data_samples$crop = substr(pair_province_prop_data_samples$sample,1,3)
  pair_province_prop_data_samples$province = str_sub(pair_province_prop_data_samples$sample,-5,-4)
  browser()
  QC_AB = pair_province_prop_data_samples %>% 
    filter(str_detect(province, "BC|AB", negate = FALSE))
  browser()
  
  group_data <- QC_AB %>%
    pair_calc_nmds() %>%
    province_treat_reps(treatment_key)
  # browser()
  province_calc_ano(QC_AB, group_data, province_title, dataset_name, outdir)
  
  return(group_data)
}

# plots the nmds for activity 2 data
plot_nmds_2 <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  
  #'*hull_var is "treatment"*
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



plot_nmds_2_all <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  
  #'*h_var is "crop"*
  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data = hull, alpha = 0.5) +
    geom_point(aes(shape = crop), size = 3) +
    aes(fill = !!hull_var) +
    labs(title = plot_title,
         x = "NMDS1", 
         y = "NMDS2",
         shape = "Crop",
         fill = tools::toTitleCase(h_var))
  
  plot
}
#'*data is province_crop_data, plot_title is NMDS: by province*
plot_nmds_2_province <- function(data, h_var, plot_title) {
  #'*Symbols are a kind of defused expression that represent objects in environments. sym() takes strings as input and turn them into symbols.*
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  # browser()
  #'*h_var is "province"*
  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2), show.legend = FALSE) + 
    geom_polygon(aes(colour = !!hull_var), data = hull, alpha = 0.5, fill = NA) +
    scale_shape_manual(values=c(7,9,10,12,3,4,8,11,0,15,1,16,2,17,25,6,5,23,60,61,35,36, 64, 38)) + #c(1, 5, 7) values=1:22
    geom_point( aes(shape = crop, colour = !!hull_var), size = 3) + #!!hull_var; show.legend = FALSE, 
    aes(fill = !!hull_var) +
    labs(title = plot_title,
         x = "NMDS1", 
         y = "NMDS2",
         shape = "Crop",
         fill = tools::toTitleCase(h_var))
  
  plot
}

#'*data is province_crop_data, plot_title is NMDS: by province*
plot_nmds_2_paired <- function(data, h_var, plot_title) {
  #'*Symbols are a kind of defused expression that represent objects in environments. sym() takes strings as input and turn them into symbols.*
  hull_var <- sym(h_var)
  
  hull <- data %>%
    group_by(!!hull_var) %>%
    slice(grDevices::chull(NMDS1, NMDS2))
  # browser()
  #'*h_var is "province"*
  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2), show.legend = FALSE) + 
    geom_polygon(aes(colour = !!hull_var), data = hull, alpha = 0.5, fill = NA) +
    scale_shape_manual(values=c(7,9,10,12,3,4,8,11,0,15,1,16,2,17,25,6,5,23,60,61,35,36, 64, 38)) + #c(1, 5, 7) values=1:22
    geom_point( aes(shape = crop, colour = !!hull_var), size = 3) + #!!hull_var; show.legend = FALSE, 
    aes(fill = !!hull_var) +
    labs(title = plot_title,
         x = "NMDS1", 
         y = "NMDS2",
         shape = "Crop",
         fill = tools::toTitleCase(h_var))
  
  plot
}

# plots and saves nmds' with separate hulls for replicate and treatment
act_1_nmds <- function(genus_data, speci_data, dataset_name, outdir) {
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
         filename = glue("{outdir}/nmds_plot_treatment_genus.png"),
         bg = "white")
  ggsave(plot = genus_reps_nmds, 
         filename = glue("{outdir}/nmds_plot_replicate_genus.png"),
         bg = "white")
  ggsave(plot = speci_treat_nmds, 
         filename = glue("{outdir}/nmds_plot_treatment_species.png"),
         bg = "white")
  ggsave(plot = speci_reps_nmds, 
         filename = glue("{outdir}/nmds_plot_replicate_species.png"),
         bg = "white")
}

# plots and saves nmds' for treatment only
# act_2_nmds <- function(genus_data, speci_data, dataset_name, outdir) {
#   genus_treat_nmds <- plot_nmds_2(genus_data,
#                                   "treatment",
#                                   "Genus NMDS - Brays Curtis")
#   speci_treat_nmds <- plot_nmds_2(speci_data,
#                                   "treatment",
#                                   "Species NMDS - Bray Curtis")
#   ggsave(plot = genus_treat_nmds,
#          filename = glue("{outdir}/nmds_plot_treatment_genus.png"),
#          bg = "white")
#   ggsave(plot = speci_treat_nmds,
#          filename = glue("{outdir}/nmds_plot_treatment_species.png"),
#          bg = "white")
# }

act_2_nmds <- function(crop_samples_data, dataset_name, outdir) {
  crop_samples_treat_nmds <- plot_nmds_2(crop_samples_data,
                                  "treatment",
                                  "Crop's Samples NMDS - Brays Curtis")
  ggsave(plot = crop_samples_treat_nmds,
         filename = glue("{outdir}/nmds_plot_treatment_crop_samples.png"),
         bg = "white")
}

all_act_2_nmds <- function(all_crop_data, dataset_name, outdir) {
  all_crops_treat_nmds <- plot_nmds_2_all(all_crop_data,
                                         "crop",
                                         "NMDS: Across Crops")
  ggsave(plot = all_crops_treat_nmds,
         filename = glue("{outdir}/nmds_plot_treatment_crop_samples.png"),
         bg = "white")
}

province_act_2_nmds <- function(province_crop_data, dataset_name, outdir) {
  province_crops_treat_nmds <- plot_nmds_2_province(province_crop_data,
                                          "province",
                                          "NMDS: by province")
  ggsave(plot = province_crops_treat_nmds,
         filename = glue("{outdir}/nmds_plot_province.png"),
         bg = "white")
}

pair_act_2_nmds <- function(province_crop_data, dataset_name, outdir) {
  province_crops_treat_nmds <- plot_nmds_2_paired(province_crop_data,
                                                    "province",
                                                    "NMDS: Paired Provinces")
  ggsave(plot = province_crops_treat_nmds,
         filename = glue("{outdir}/nmds_plot_province.png"),
         bg = "white")
}


# make and save nmds plots for genus and species levels using read data
# uses raw reads to calculate proportions
export("make_nmds_plots")
make_nmds_plots <- function(data, treatment_key, dataset_name, outdir) {
  file.create(glue("{outdir}/anosim.txt"))
  # genus_data <- filter(data, taxRank == "G") %>%
  #   prep_and_ano(treatment_key, "Genus", dataset_name, outdir)
  # speci_data <- filter(data, taxRank == "S") %>%
  #   prep_and_ano(treatment_key, "Species", dataset_name, outdir)
  
  crop_samples_data <- prep_and_ano(data, treatment_key, "Crop's Samples", dataset_name, outdir)
  # browser()
  # crop_samples_data <- select(crop_samples_data, -"crop")
  # browser()
  utils::write.csv(crop_samples_data,
                   file = glue("{outdir}/nmds_plot_data_crop_samples.csv"),
                   row.names = F)
  # utils::write.csv(speci_data,
  #                  file = glue("{outdir}/nmds_plot_data_species.csv"),
  #                  row.names = F)
  
  # split act 1 and 2 nmds here, based on number of treatments
  # if (length(unique(genus_data$treatment)) > 2) {
  #   # can make convex hulls for both replicate and treatment
  #   act_1_nmds(genus_data, speci_data, dataset_name, outdir)
  # } else {
  #   # can only make hulls for treatment
  #   act_2_nmds(genus_data, speci_data, dataset_name, outdir)
  # }

    # can only make hulls for treatment
  act_2_nmds(crop_samples_data, dataset_name, outdir)
}

export("make_nmds_plots_all_crops")
make_nmds_plots_all_crops <- function(merged_crops, treatment_key, dataset_name, outdir) {
  file.create(glue("{outdir}/anosim.txt"))
 
  all_crop_data <- all_prep_and_ano(merged_crops, treatment_key, "All Crops Compared", dataset_name, outdir)
  #browser()
  
  utils::write.csv(all_crop_data,
                   file = glue("{outdir}/nmds_plot_data_all_crops.csv"),
                   row.names = F)
 
  all_act_2_nmds(all_crop_data, dataset_name, outdir)
}

export("make_nmds_plots_provinces")
make_nmds_plots_provinces <- function(merged_crops, treatment_key, dataset_name, outdir) {
  browser()
  #'*at this point HBB21 is already corrupted - in the form HBB01u_t2>, when all the other provinces are in the form COR04e_t2_ON_20*
  province_data <- province_prep_and_ano(merged_crops, treatment_key, "provinces Compared", dataset_name, outdir)
  
  
  utils::write.csv(province_data,
                   file = glue("{outdir}/nmds_plot_data_provinces.csv"),
                   row.names = F)
  
  province_act_2_nmds(province_data, dataset_name, outdir)
}

export("make_nmds_plots_paired_provinces")
make_nmds_plots_paired_provinces <- function(merged_crops, treatment_key, dataset_name, outdir) {
  # browser()
  #'*at this point HBB21 is already corrupted - in the form HBB01u_t2>, when all the other provinces are in the form COR04e_t2_ON_20*
  BC_AB_paired <- BC_AB_pair_prep_and_ano(merged_crops, treatment_key, "provinces Compared", dataset_name, outdir)
  
  utils::write.csv(BC_AB_paired,
                   file = glue("{outdir}/nmds_plot_data_BC_AB.csv"),
                   row.names = F)
  
  pair_act_2_nmds(BC_AB_paired, dataset_name, outdir)
  
  QC_AB_paired <- QC_AB_pair_prep_and_ano(merged_crops, treatment_key, "provinces Compared", dataset_name, outdir)
  
  utils::write.csv(QC_AB_paired,
                   file = glue("{outdir}/nmds_plot_data_QC_AB.csv"),
                   row.names = F)
  
  pair_act_2_nmds(QC_AB_paired, dataset_name, outdir)
}