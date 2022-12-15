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
# interest_list <- c("Tetracyclines",
#                    "Sulfonamides",
#                    "Rifampin",
#                    "Mupirocin",
#                    "Multi-metal_resistance",
#                    "Multi-drug_resistance",
#                    "Multi-biocide_resistance",
#                    "MLS",
#                    "Iron_resistance",
#                    "Glycopeptides",
#                    "Drug_and_biocide_resistance",
#                    "Drug_and_biocide_and_metal_resistance",
#                    "Copper_resistance",
#                    "betalactams",
#                    "Aminoglycosides",
#                    "Aminocoumarins")


# Wrangling Functions -----------------------------------------------------
# cleans data into tidy format
#'*commented out to preserve original code*'
#'#'*select() is removing those coluMBs, which we don't have in the first place*
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
  crop_label = substr(clean_data$sample,1,3)
  year_label = str_sub(clean_data$sample, -2, -1)
  #browser()
  clean_data$crop = paste(crop_label,year_label, sep = "")
  # browser()

    # data.frame(crop = c("HBB"))
  clean_data = ungroup(clean_data)
    return(clean_data)
}

tidy_data_RA_samples <- function(data) {
  clean_data <- select(data, -geneID, -type, -mechanism, -gene) %>%
    # browser()
    pivot_longer(!class, names_to = "sample", values_to = "value")
  # browser()
  clean_data[is.na(clean_data)] <- 0
  clean_data = group_by(clean_data, sample, class) %>%
    summarise(value=sum(value)) %>%
    # browser()
    pivot_wider(names_from = "class", values_from = "value")
  crop_label = substr(clean_data$sample,1,3)
  year_label = str_sub(clean_data$sample, -2, -1)
  # browser()
  clean_data$crop = paste(crop_label,year_label, sep = "")
  # browser()
  
  # data.frame(crop = c("HBB"))
  clean_data = ungroup(clean_data)
  colnames(clean_data)[37] = "Quaternary_Ammonium_Compounds_resistance"
  colnames(clean_data)[26] = "Multi_biocide_resistance"
  colnames(clean_data)[27] = "Multi_drug_resistance"
  colnames(clean_data)[28] = "Multi_metal_resistance"
  browser()
  
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
  # browser()
  crop_label = substr(clean_data$sample,1,3)
  year_label = str_sub(clean_data$sample, -2, -1)
  treatment_label = substr(clean_data$sample,6,6)
  # browser()
  clean_data$crop = paste(crop_label,year_label,treatment_label, sep = "")
  clean_data$province = str_sub(clean_data$sample, -5, -4)
  # browser()
  
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

#'*use for slimming down ARGs in legends*
# calculates and adds Other Bacteria col
#'*d is plot_data <- filter(data, name %in% interest_list)*
add_other_bac <- function(d) {
  coluMB_index <- !(names(d) %in% "class")
  browser()
  # Get the sum counts for taxa of interest, except Bacteria (domain), then
  # subtract their sum from Bacteria clade count. This leaves the Bacteria row
  # representing all taxa of non-interest
  df <- d[d$name!="Bacteria", coluMB_index]
  d[d$name=="Bacteria", coluMB_index] <-
    d[d$name=="Bacteria", coluMB_index] - as.list(colSums(df, na.rm=T))
  # Change the grouping name
  d[d$name=="Bacteria",]$name <- "Other Bacteria"

  return(d)
}

# calculates and adds Other Bacteria col
# add_other_bac <- function(d) {
#   coluMB_index <- !(names(d) %in% c("taxRank","taxLineage","taxID","depth","name"))
#   # Get the sum counts for taxa of interest, except Bacteria (domain), then
#   # subtract their sum from Bacteria clade count. This leaves the Bacteria row
#   # representing all taxa of non-interest
#   df <- d[d$name!="Bacteria", coluMB_index]
#   d[d$name=="Bacteria", coluMB_index] <-
#     d[d$name=="Bacteria", coluMB_index] - as.list(colSums(df, na.rm=T))
#   # Change the grouping name
#   d[d$name=="Bacteria",]$name <- "Other ARGs"
# 
#   return(d)
# }



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
  crop_col <- select(d, "crop")
  plot_data_samples <- select(d, -"sample", -"crop") %>%
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
  # browser()
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample) %>%
    mutate(crop_col)


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
  # browser()
  sample_col <- select(d, "sample")
  crop_col <- select(d, "crop")
  province_col = select(d, "province")
  # browser()
  plot_data_samples <- select(d, -"sample", -"crop", -"province") %>%
    apply(MARGIN = 1,
          FUN = function(x) x / sum(x) * 100) %>%
    #browser()
    t() %>%
    as.data.frame() %>%
    mutate(sample_col) %>%
    relocate(sample) %>% 
    mutate(crop_col) %>%
    relocate(crop) %>% 
    mutate(province_col) %>%
      relocate(province)

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
# Samples gathered from coluMB names are treated with regular expressions to
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
    stop("The sample coluMB names do not match the implemented regex patterns")
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
                            cols = !c("sample", "treatment", "replicate", "crop", "province"), #"crop"
                            names_to = "ARGs",
                            values_to = "value")
  # browser()
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

# Reorders table using taxa coluMB
order_taxa <- function(d) {
  # browser()
  taxa_order <- get_taxa_order(d)
  d$ARGs <- factor(d$ARGs, levels=taxa_order$Group.1)
  return(d)
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
  # browser()
  # d$crop = substr(d$sample,1,3)
  # browser()
  pd_avg <- aggregate(d[,c("value")], list(d$ARGs), mean) %>%
    arrange(value)
  # browser()
  #taxa_order <- append(pd_avg$Group.1, "Other Bacteria")
}


# Relative Abundance ------------------------------------------------------

# plots relative abundance data for taxa of interest
# uses the taxa, value, treatment, and replicate coluMB from data
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
                           aes(x = crop,
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
  # browser()
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
export("prov_make_interest_abundance")
prov_make_interest_abundance <- function(data, treatment_key, dataset_name, outdir) {  #additional_taxa,

  plot_data_samples <- tidy_data_province(data) 
  # browser()
  
  AB_RA_data = plot_data_samples %>% 
    filter(str_detect(province, "AB"))
  # browser()
  AB_RA_data$Acetate_resistance = sum(AB_RA_data$Acetate_resistance)
  AB_RA_data$Acid_resistance = sum(AB_RA_data$Acid_resistance)
  AB_RA_data$Aldehyde_resistance = sum(AB_RA_data$Aldehyde_resistance)
  AB_RA_data$Aluminum_resistance = sum(AB_RA_data$Aluminum_resistance)
  AB_RA_data$Aminocoumarins = sum(AB_RA_data$Aminocoumarins)
  AB_RA_data$Aminoglycosides = sum(AB_RA_data$Aminoglycosides)
  AB_RA_data$Bacitracin = sum(AB_RA_data$Bacitracin)
  AB_RA_data$Biguanide_resistance = sum(AB_RA_data$Biguanide_resistance)
  AB_RA_data$betalactams = sum(AB_RA_data$betalactams)
  AB_RA_data$Cadmium_resistance = sum(AB_RA_data$Cadmium_resistance)
  AB_RA_data$Cationic_antimicrobial_peptides = sum(AB_RA_data$Cationic_antimicrobial_peptides)
  AB_RA_data$Chromium_resistance = sum(AB_RA_data$Chromium_resistance)
  AB_RA_data$Cobalt_resistance = sum(AB_RA_data$Cobalt_resistance)
  AB_RA_data$Copper_resistance = sum(AB_RA_data$Copper_resistance)

  AB_RA_data$Drug_and_biocide_and_metal_resistance = sum(AB_RA_data$Drug_and_biocide_and_metal_resistance)
  AB_RA_data$Drug_and_biocide_resistance = sum(AB_RA_data$Drug_and_biocide_resistance)
  AB_RA_data$Drug_and_metal_resistance = sum(AB_RA_data$Drug_and_metal_resistance)
  AB_RA_data$Fluoroquinolones = sum(AB_RA_data$Fluoroquinolones)
  AB_RA_data$Fosfomycin = sum(AB_RA_data$Fosfomycin)
  AB_RA_data$Fusidic_acid = sum(AB_RA_data$Fusidic_acid)
  AB_RA_data$Glycopeptides = sum(AB_RA_data$Glycopeptides)
  AB_RA_data$Iron_resistance = sum(AB_RA_data$Iron_resistance)
  AB_RA_data$Lipopeptides = sum(AB_RA_data$Lipopeptides)
  AB_RA_data$Mercury_resistance = sum(AB_RA_data$Mercury_resistance)
  AB_RA_data$Metronidazole = sum(AB_RA_data$Metronidazole)
  AB_RA_data$MLS = sum(AB_RA_data$MLS)
  AB_RA_data$Multi_biocide_resistance = sum(AB_RA_data$Multi_biocide_resistance)
  AB_RA_data$Multi_drug_resistance = sum(AB_RA_data$Multi_drug_resistance)
  
  AB_RA_data$Multi_metal_resistance = sum(AB_RA_data$Multi_metal_resistance)
  AB_RA_data$Mupirocin = sum(AB_RA_data$Mupirocin)
  AB_RA_data$Naphthoquinone = sum(AB_RA_data$Naphthoquinone)
  AB_RA_data$Mycobacterium_tuberculosis_specific_Drug = sum(AB_RA_data$Mycobacterium_tuberculosis_specific_Drug)
  AB_RA_data$Nickel_resistance = sum(AB_RA_data$Nickel_resistance)
  AB_RA_data$Nucleosides = sum(AB_RA_data$Nucleosides)
  AB_RA_data$Paraquat_resistance = sum(AB_RA_data$Paraquat_resistance)
  AB_RA_data$Peroxide_resistance = sum(AB_RA_data$Peroxide_resistance)
  AB_RA_data$Phenicol = sum(AB_RA_data$Phenicol)
  AB_RA_data$Phenolic_compound_resistance = sum(AB_RA_data$Phenolic_compound_resistance)
  AB_RA_data$Polyamine_resistance = sum(AB_RA_data$Polyaqacmine_resistance)
  AB_RA_data$QACs_resistance = sum(AB_RA_data$QACs_resistance)
  AB_RA_data$Rifampin = sum(AB_RA_data$Rifampin)
  AB_RA_data$Sodium_resistance = sum(AB_RA_data$Sodium_resistance)
  
  AB_RA_data$Sulfonamides = sum(AB_RA_data$Sulfonamides)
  AB_RA_data$Tellurium_resistance = sum(AB_RA_data$Tellurium_resistance)
  AB_RA_data$Tetracenomycin = sum(AB_RA_data$Tetracenomycin)
  AB_RA_data$Tetracyclines = sum(AB_RA_data$Tetracyclines)
  AB_RA_data$Trimethoprim = sum(AB_RA_data$Trimethoprim)
  AB_RA_data$Tungsten_Resistance = sum(AB_RA_data$Tungsten_Resistance)
  AB_RA_data$Zinc_resistance = sum(AB_RA_data$Zinc_resistance)
  
  AB_RA_data = AB_RA_data[1:1,]
  # other_ARGs select(AB_RA_data, !c("Tetracyclines"))

  # browser()
  plot_data_AB_RA = calc_prop_province(AB_RA_data) %>%
    treat_reps(treatment_key)
  # browser()
  #' #'*sends the relative abundance data up to the ggplot function*
  plot_AB_RA <- tidy_to_long_samples(plot_data_AB_RA)
  # browser()
  plot_AB_RA <- order_taxa(plot_AB_RA)
  # browser()
  plot_AB_RA <- plot_interest_abundance_crop_av(plot_AB_RA)

  ggsave(plot = plot_AB_RA,
         filename = glue("{outdir}/relative_abundance_AB.png"),
         bg = "white")
  
  utils::write.csv(plot_data_AB_RA,
                   file = glue("{outdir}/relative_proportions_abundance_AB.csv"),
                   row.names = F)

#######
  ON_RA_data = plot_data_samples %>% 
    filter(str_detect(province, "ON"))
  # browser()
  ON_RA_data$Acetate_resistance = sum(ON_RA_data$Acetate_resistance)
  ON_RA_data$Acid_resistance = sum(ON_RA_data$Acid_resistance)
  ON_RA_data$Aldehyde_resistance = sum(ON_RA_data$Aldehyde_resistance)
  ON_RA_data$Aluminum_resistance = sum(ON_RA_data$Aluminum_resistance)
  ON_RA_data$Aminocoumarins = sum(ON_RA_data$Aminocoumarins)
  ON_RA_data$Aminoglycosides = sum(ON_RA_data$Aminoglycosides)
  ON_RA_data$Bacitracin = sum(ON_RA_data$Bacitracin)
  ON_RA_data$Biguanide_resistance = sum(ON_RA_data$Biguanide_resistance)
  ON_RA_data$betalactams = sum(ON_RA_data$betalactams)
  ON_RA_data$Cadmium_resistance = sum(ON_RA_data$Cadmium_resistance)
  ON_RA_data$Cationic_antimicrobial_peptides = sum(ON_RA_data$Cationic_antimicrobial_peptides)
  ON_RA_data$Chromium_resistance = sum(ON_RA_data$Chromium_resistance)
  ON_RA_data$Cobalt_resistance = sum(ON_RA_data$Cobalt_resistance)
  ON_RA_data$Copper_resistance = sum(ON_RA_data$Copper_resistance)
  
  ON_RA_data$Drug_and_biocide_and_metal_resistance = sum(ON_RA_data$Drug_and_biocide_and_metal_resistance)
  ON_RA_data$Drug_and_biocide_resistance = sum(ON_RA_data$Drug_and_biocide_resistance)
  ON_RA_data$Drug_and_metal_resistance = sum(ON_RA_data$Drug_and_metal_resistance)
  ON_RA_data$Fluoroquinolones = sum(ON_RA_data$Fluoroquinolones)
  ON_RA_data$Fosfomycin = sum(ON_RA_data$Fosfomycin)
  ON_RA_data$Fusidic_acid = sum(ON_RA_data$Fusidic_acid)
  ON_RA_data$Glycopeptides = sum(ON_RA_data$Glycopeptides)
  ON_RA_data$Iron_resistance = sum(ON_RA_data$Iron_resistance)
  ON_RA_data$Lipopeptides = sum(ON_RA_data$Lipopeptides)
  ON_RA_data$Mercury_resistance = sum(ON_RA_data$Mercury_resistance)
  ON_RA_data$Metronidazole = sum(ON_RA_data$Metronidazole)
  ON_RA_data$MLS = sum(ON_RA_data$MLS)
  ON_RA_data$Multi_biocide_resistance = sum(ON_RA_data$Multi_biocide_resistance)
  ON_RA_data$Multi_drug_resistance = sum(ON_RA_data$Multi_drug_resistance)
  
  ON_RA_data$Multi_metal_resistance = sum(ON_RA_data$Multi_metal_resistance)
  ON_RA_data$Mupirocin = sum(ON_RA_data$Mupirocin)
  ON_RA_data$Naphthoquinone = sum(ON_RA_data$Naphthoquinone)
  ON_RA_data$Mycobacterium_tuberculosis_specific_Drug = sum(ON_RA_data$Mycobacterium_tuberculosis_specific_Drug)
  ON_RA_data$Nickel_resistance = sum(ON_RA_data$Nickel_resistance)
  ON_RA_data$Nucleosides = sum(ON_RA_data$Nucleosides)
  ON_RA_data$Paraquat_resistance = sum(ON_RA_data$Paraquat_resistance)
  ON_RA_data$Peroxide_resistance = sum(ON_RA_data$Peroxide_resistance)
  ON_RA_data$Phenicol = sum(ON_RA_data$Phenicol)
  ON_RA_data$Phenolic_compound_resistance = sum(ON_RA_data$Phenolic_compound_resistance)
  ON_RA_data$Polyamine_resistance = sum(ON_RA_data$Polyaqacmine_resistance)
  ON_RA_data$QACs_resistance = sum(ON_RA_data$QACs_resistance)
  ON_RA_data$Rifampin = sum(ON_RA_data$Rifampin)
  ON_RA_data$Sodium_resistance = sum(ON_RA_data$Sodium_resistance)
  
  ON_RA_data$Sulfonamides = sum(ON_RA_data$Sulfonamides)
  ON_RA_data$Tellurium_resistance = sum(ON_RA_data$Tellurium_resistance)
  ON_RA_data$Tetracenomycin = sum(ON_RA_data$Tetracenomycin)
  ON_RA_data$Tetracyclines = sum(ON_RA_data$Tetracyclines)
  ON_RA_data$Trimethoprim = sum(ON_RA_data$Trimethoprim)
  ON_RA_data$Tungsten_Resistance = sum(ON_RA_data$Tungsten_Resistance)
  ON_RA_data$Zinc_resistance = sum(ON_RA_data$Zinc_resistance)
  
  ON_RA_data = ON_RA_data[1:1,]
  # other_ARGs select(ON_RA_data, !c("Tetracyclines"))
  
  # browser()
  plot_data_ON_RA = calc_prop_province(ON_RA_data) %>%
    treat_reps(treatment_key)
  # browser()
  #' #'*sends the relative abundance data up to the ggplot function*
  plot_ON_RA <- tidy_to_long_samples(plot_data_ON_RA)
  # browser()
  plot_ON_RA <- order_taxa(plot_ON_RA)
  # browser()
  plot_ON_RA <- plot_interest_abundance_crop_av(plot_ON_RA)
  
  ggsave(plot = plot_ON_RA,
         filename = glue("{outdir}/relative_abundance_ON.png"),
         bg = "white")
  
  utils::write.csv(plot_data_ON_RA,
                   file = glue("{outdir}/relative_proportions_abundance_ON.csv"),
                   row.names = F)
  
  #########
  
  MB_RA_data = plot_data_samples %>% 
    filter(str_detect(province, "MB"))
  # browser()
  MB_RA_data$Acetate_resistance = sum(MB_RA_data$Acetate_resistance)
  MB_RA_data$Acid_resistance = sum(MB_RA_data$Acid_resistance)
  MB_RA_data$Aldehyde_resistance = sum(MB_RA_data$Aldehyde_resistance)
  MB_RA_data$Aluminum_resistance = sum(MB_RA_data$Aluminum_resistance)
  MB_RA_data$Aminocoumarins = sum(MB_RA_data$Aminocoumarins)
  MB_RA_data$Aminoglycosides = sum(MB_RA_data$Aminoglycosides)
  MB_RA_data$Bacitracin = sum(MB_RA_data$Bacitracin)
  MB_RA_data$Biguanide_resistance = sum(MB_RA_data$Biguanide_resistance)
  MB_RA_data$betalactams = sum(MB_RA_data$betalactams)
  MB_RA_data$Cadmium_resistance = sum(MB_RA_data$Cadmium_resistance)
  MB_RA_data$Cationic_antimicrobial_peptides = sum(MB_RA_data$Cationic_antimicrobial_peptides)
  MB_RA_data$Chromium_resistance = sum(MB_RA_data$Chromium_resistance)
  MB_RA_data$Cobalt_resistance = sum(MB_RA_data$Cobalt_resistance)
  MB_RA_data$Copper_resistance = sum(MB_RA_data$Copper_resistance)
  
  MB_RA_data$Drug_and_biocide_and_metal_resistance = sum(MB_RA_data$Drug_and_biocide_and_metal_resistance)
  MB_RA_data$Drug_and_biocide_resistance = sum(MB_RA_data$Drug_and_biocide_resistance)
  MB_RA_data$Drug_and_metal_resistance = sum(MB_RA_data$Drug_and_metal_resistance)
  MB_RA_data$Fluoroquinolones = sum(MB_RA_data$Fluoroquinolones)
  MB_RA_data$Fosfomycin = sum(MB_RA_data$Fosfomycin)
  MB_RA_data$Fusidic_acid = sum(MB_RA_data$Fusidic_acid)
  MB_RA_data$Glycopeptides = sum(MB_RA_data$Glycopeptides)
  MB_RA_data$Iron_resistance = sum(MB_RA_data$Iron_resistance)
  MB_RA_data$Lipopeptides = sum(MB_RA_data$Lipopeptides)
  MB_RA_data$Mercury_resistance = sum(MB_RA_data$Mercury_resistance)
  MB_RA_data$Metronidazole = sum(MB_RA_data$Metronidazole)
  MB_RA_data$MLS = sum(MB_RA_data$MLS)
  MB_RA_data$Multi_biocide_resistance = sum(MB_RA_data$Multi_biocide_resistance)
  MB_RA_data$Multi_drug_resistance = sum(MB_RA_data$Multi_drug_resistance)
  
  MB_RA_data$Multi_metal_resistance = sum(MB_RA_data$Multi_metal_resistance)
  MB_RA_data$Mupirocin = sum(MB_RA_data$Mupirocin)
  MB_RA_data$Naphthoquinone = sum(MB_RA_data$Naphthoquinone)
  MB_RA_data$Mycobacterium_tuberculosis_specific_Drug = sum(MB_RA_data$Mycobacterium_tuberculosis_specific_Drug)
  MB_RA_data$Nickel_resistance = sum(MB_RA_data$Nickel_resistance)
  MB_RA_data$Nucleosides = sum(MB_RA_data$Nucleosides)
  MB_RA_data$Paraquat_resistance = sum(MB_RA_data$Paraquat_resistance)
  MB_RA_data$Peroxide_resistance = sum(MB_RA_data$Peroxide_resistance)
  MB_RA_data$Phenicol = sum(MB_RA_data$Phenicol)
  MB_RA_data$Phenolic_compound_resistance = sum(MB_RA_data$Phenolic_compound_resistance)
  MB_RA_data$Polyamine_resistance = sum(MB_RA_data$Polyaqacmine_resistance)
  MB_RA_data$QACs_resistance = sum(MB_RA_data$QACs_resistance)
  MB_RA_data$Rifampin = sum(MB_RA_data$Rifampin)
  MB_RA_data$Sodium_resistance = sum(MB_RA_data$Sodium_resistance)
  
  MB_RA_data$Sulfonamides = sum(MB_RA_data$Sulfonamides)
  MB_RA_data$Tellurium_resistance = sum(MB_RA_data$Tellurium_resistance)
  MB_RA_data$Tetracenomycin = sum(MB_RA_data$Tetracenomycin)
  MB_RA_data$Tetracyclines = sum(MB_RA_data$Tetracyclines)
  MB_RA_data$Trimethoprim = sum(MB_RA_data$Trimethoprim)
  MB_RA_data$Tungsten_Resistance = sum(MB_RA_data$Tungsten_Resistance)
  MB_RA_data$Zinc_resistance = sum(MB_RA_data$Zinc_resistance)
  
  MB_RA_data = MB_RA_data[1:1,]
  # other_ARGs select(MB_RA_data, !c("Tetracyclines"))
  
  # browser()
  plot_data_MB_RA = calc_prop_province(MB_RA_data) %>%
    treat_reps(treatment_key)
  # browser()
  #' #'*sends the relative abundance data up to the ggplot function*
  plot_MB_RA <- tidy_to_long_samples(plot_data_MB_RA)
  # browser()
  plot_MB_RA <- order_taxa(plot_MB_RA)
  # browser()
  plot_MB_RA <- plot_interest_abundance_crop_av(plot_MB_RA)
  
  ggsave(plot = plot_MB_RA,
         filename = glue("{outdir}/relative_abundance_MB.png"),
         bg = "white")
  
  utils::write.csv(plot_data_MB_RA,
                   file = glue("{outdir}/relative_proportions_abundance_MB.csv"),
                   row.names = F)
  #########
  BC_RA_data = plot_data_samples %>% 
    filter(str_detect(province, "BC"))
  # browser()
  BC_RA_data$Acetate_resistance = sum(BC_RA_data$Acetate_resistance)
  BC_RA_data$Acid_resistance = sum(BC_RA_data$Acid_resistance)
  BC_RA_data$Aldehyde_resistance = sum(BC_RA_data$Aldehyde_resistance)
  BC_RA_data$Aluminum_resistance = sum(BC_RA_data$Aluminum_resistance)
  BC_RA_data$Aminocoumarins = sum(BC_RA_data$Aminocoumarins)
  BC_RA_data$Aminoglycosides = sum(BC_RA_data$Aminoglycosides)
  BC_RA_data$Bacitracin = sum(BC_RA_data$Bacitracin)
  BC_RA_data$Biguanide_resistance = sum(BC_RA_data$Biguanide_resistance)
  BC_RA_data$betalactams = sum(BC_RA_data$betalactams)
  BC_RA_data$Cadmium_resistance = sum(BC_RA_data$Cadmium_resistance)
  BC_RA_data$Cationic_antimicrobial_peptides = sum(BC_RA_data$Cationic_antimicrobial_peptides)
  BC_RA_data$Chromium_resistance = sum(BC_RA_data$Chromium_resistance)
  BC_RA_data$Cobalt_resistance = sum(BC_RA_data$Cobalt_resistance)
  BC_RA_data$Copper_resistance = sum(BC_RA_data$Copper_resistance)
  
  BC_RA_data$Drug_and_biocide_and_metal_resistance = sum(BC_RA_data$Drug_and_biocide_and_metal_resistance)
  BC_RA_data$Drug_and_biocide_resistance = sum(BC_RA_data$Drug_and_biocide_resistance)
  BC_RA_data$Drug_and_metal_resistance = sum(BC_RA_data$Drug_and_metal_resistance)
  BC_RA_data$Fluoroquinolones = sum(BC_RA_data$Fluoroquinolones)
  BC_RA_data$Fosfomycin = sum(BC_RA_data$Fosfomycin)
  BC_RA_data$Fusidic_acid = sum(BC_RA_data$Fusidic_acid)
  BC_RA_data$Glycopeptides = sum(BC_RA_data$Glycopeptides)
  BC_RA_data$Iron_resistance = sum(BC_RA_data$Iron_resistance)
  BC_RA_data$Lipopeptides = sum(BC_RA_data$Lipopeptides)
  BC_RA_data$Mercury_resistance = sum(BC_RA_data$Mercury_resistance)
  BC_RA_data$Metronidazole = sum(BC_RA_data$Metronidazole)
  BC_RA_data$MLS = sum(BC_RA_data$MLS)
  BC_RA_data$Multi_biocide_resistance = sum(BC_RA_data$Multi_biocide_resistance)
  BC_RA_data$Multi_drug_resistance = sum(BC_RA_data$Multi_drug_resistance)
  
  BC_RA_data$Multi_metal_resistance = sum(BC_RA_data$Multi_metal_resistance)
  BC_RA_data$Mupirocin = sum(BC_RA_data$Mupirocin)
  BC_RA_data$Naphthoquinone = sum(BC_RA_data$Naphthoquinone)
  BC_RA_data$Mycobacterium_tuberculosis_specific_Drug = sum(BC_RA_data$Mycobacterium_tuberculosis_specific_Drug)
  BC_RA_data$Nickel_resistance = sum(BC_RA_data$Nickel_resistance)
  BC_RA_data$Nucleosides = sum(BC_RA_data$Nucleosides)
  BC_RA_data$Paraquat_resistance = sum(BC_RA_data$Paraquat_resistance)
  BC_RA_data$Peroxide_resistance = sum(BC_RA_data$Peroxide_resistance)
  BC_RA_data$Phenicol = sum(BC_RA_data$Phenicol)
  BC_RA_data$Phenolic_compound_resistance = sum(BC_RA_data$Phenolic_compound_resistance)
  BC_RA_data$Polyamine_resistance = sum(BC_RA_data$Polyaqacmine_resistance)
  BC_RA_data$QACs_resistance = sum(BC_RA_data$QACs_resistance)
  BC_RA_data$Rifampin = sum(BC_RA_data$Rifampin)
  BC_RA_data$Sodium_resistance = sum(BC_RA_data$Sodium_resistance)
  
  BC_RA_data$Sulfonamides = sum(BC_RA_data$Sulfonamides)
  BC_RA_data$Tellurium_resistance = sum(BC_RA_data$Tellurium_resistance)
  BC_RA_data$Tetracenomycin = sum(BC_RA_data$Tetracenomycin)
  BC_RA_data$Tetracyclines = sum(BC_RA_data$Tetracyclines)
  BC_RA_data$Trimethoprim = sum(BC_RA_data$Trimethoprim)
  BC_RA_data$Tungsten_Resistance = sum(BC_RA_data$Tungsten_Resistance)
  BC_RA_data$Zinc_resistance = sum(BC_RA_data$Zinc_resistance)
  
  BC_RA_data = BC_RA_data[1:1,]
  # other_ARGs select(BC_RA_data, !c("Tetracyclines"))
  
  # browser()
  plot_data_BC_RA = calc_prop_province(BC_RA_data) %>%
    treat_reps(treatment_key)
  # browser()
  #' #'*sends the relative abundance data up to the ggplot function*
  plot_BC_RA <- tidy_to_long_samples(plot_data_BC_RA)
  # browser()
  plot_BC_RA <- order_taxa(plot_BC_RA)
  # browser()
  plot_BC_RA <- plot_interest_abundance_crop_av(plot_BC_RA)
  
  ggsave(plot = plot_BC_RA,
         filename = glue("{outdir}/relative_abundance_BC.png"),
         bg = "white")
  
  utils::write.csv(plot_data_BC_RA,
                   file = glue("{outdir}/relative_proportions_abundance_BC.csv"),
                   row.names = F)
  
  #######
  QB_RA_data = plot_data_samples %>% 
    filter(str_detect(province, "QB"))
  # browser()
  QB_RA_data$Acetate_resistance = sum(QB_RA_data$Acetate_resistance)
  QB_RA_data$Acid_resistance = sum(QB_RA_data$Acid_resistance)
  QB_RA_data$Aldehyde_resistance = sum(QB_RA_data$Aldehyde_resistance)
  QB_RA_data$Aluminum_resistance = sum(QB_RA_data$Aluminum_resistance)
  QB_RA_data$Aminocoumarins = sum(QB_RA_data$Aminocoumarins)
  QB_RA_data$Aminoglycosides = sum(QB_RA_data$Aminoglycosides)
  QB_RA_data$Bacitracin = sum(QB_RA_data$Bacitracin)
  QB_RA_data$Biguanide_resistance = sum(QB_RA_data$Biguanide_resistance)
  QB_RA_data$betalactams = sum(QB_RA_data$betalactams)
  QB_RA_data$Cadmium_resistance = sum(QB_RA_data$Cadmium_resistance)
  QB_RA_data$Cationic_antimicrobial_peptides = sum(QB_RA_data$Cationic_antimicrobial_peptides)
  QB_RA_data$Chromium_resistance = sum(QB_RA_data$Chromium_resistance)
  QB_RA_data$Cobalt_resistance = sum(QB_RA_data$Cobalt_resistance)
  QB_RA_data$Copper_resistance = sum(QB_RA_data$Copper_resistance)
  
  QB_RA_data$Drug_and_biocide_and_metal_resistance = sum(QB_RA_data$Drug_and_biocide_and_metal_resistance)
  QB_RA_data$Drug_and_biocide_resistance = sum(QB_RA_data$Drug_and_biocide_resistance)
  QB_RA_data$Drug_and_metal_resistance = sum(QB_RA_data$Drug_and_metal_resistance)
  QB_RA_data$Fluoroquinolones = sum(QB_RA_data$Fluoroquinolones)
  QB_RA_data$Fosfomycin = sum(QB_RA_data$Fosfomycin)
  QB_RA_data$Fusidic_acid = sum(QB_RA_data$Fusidic_acid)
  QB_RA_data$Glycopeptides = sum(QB_RA_data$Glycopeptides)
  QB_RA_data$Iron_resistance = sum(QB_RA_data$Iron_resistance)
  QB_RA_data$Lipopeptides = sum(QB_RA_data$Lipopeptides)
  QB_RA_data$Mercury_resistance = sum(QB_RA_data$Mercury_resistance)
  QB_RA_data$Metronidazole = sum(QB_RA_data$Metronidazole)
  QB_RA_data$MLS = sum(QB_RA_data$MLS)
  QB_RA_data$Multi_biocide_resistance = sum(QB_RA_data$Multi_biocide_resistance)
  QB_RA_data$Multi_drug_resistance = sum(QB_RA_data$Multi_drug_resistance)
  
  QB_RA_data$Multi_metal_resistance = sum(QB_RA_data$Multi_metal_resistance)
  QB_RA_data$Mupirocin = sum(QB_RA_data$Mupirocin)
  QB_RA_data$Naphthoquinone = sum(QB_RA_data$Naphthoquinone)
  QB_RA_data$Mycobacterium_tuberculosis_specific_Drug = sum(QB_RA_data$Mycobacterium_tuberculosis_specific_Drug)
  QB_RA_data$Nickel_resistance = sum(QB_RA_data$Nickel_resistance)
  QB_RA_data$Nucleosides = sum(QB_RA_data$Nucleosides)
  QB_RA_data$Paraquat_resistance = sum(QB_RA_data$Paraquat_resistance)
  QB_RA_data$Peroxide_resistance = sum(QB_RA_data$Peroxide_resistance)
  QB_RA_data$Phenicol = sum(QB_RA_data$Phenicol)
  QB_RA_data$Phenolic_compound_resistance = sum(QB_RA_data$Phenolic_compound_resistance)
  QB_RA_data$Polyamine_resistance = sum(QB_RA_data$Polyaqacmine_resistance)
  QB_RA_data$QACs_resistance = sum(QB_RA_data$QACs_resistance)
  QB_RA_data$Rifampin = sum(QB_RA_data$Rifampin)
  QB_RA_data$Sodium_resistance = sum(QB_RA_data$Sodium_resistance)
  
  QB_RA_data$Sulfonamides = sum(QB_RA_data$Sulfonamides)
  QB_RA_data$Tellurium_resistance = sum(QB_RA_data$Tellurium_resistance)
  QB_RA_data$Tetracenomycin = sum(QB_RA_data$Tetracenomycin)
  QB_RA_data$Tetracyclines = sum(QB_RA_data$Tetracyclines)
  QB_RA_data$Trimethoprim = sum(QB_RA_data$Trimethoprim)
  QB_RA_data$Tungsten_Resistance = sum(QB_RA_data$Tungsten_Resistance)
  QB_RA_data$Zinc_resistance = sum(QB_RA_data$Zinc_resistance)
  
  QB_RA_data = QB_RA_data[1:1,]
  # other_ARGs select(QB_RA_data, !c("Tetracyclines"))
  
  # browser()
  plot_data_QB_RA = calc_prop_province(QB_RA_data) %>%
    treat_reps(treatment_key)
  browser()
  #' #'*sends the relative abundance data up to the ggplot function*
  plot_QB_RA <- tidy_to_long_samples(plot_data_QB_RA)
  # browser()
  plot_QB_RA <- order_taxa(plot_QB_RA)
  # browser()
  plot_QB_RA <- plot_interest_abundance_crop_av(plot_QB_RA)
  
  ggsave(plot = plot_QB_RA,
         filename = glue("{outdir}/relative_abundance_QB.png"),
         bg = "white")
  
  utils::write.csv(plot_data_QB_RA,
                   file = glue("{outdir}/relative_proportions_abundance_QB.csv"),
                   row.names = F)
  
}

# Returns relative abundance for taxa of interest
#'*starting point, that calls all other functions in this script*
#'*"data" is being fed "tables[["raw_clade"]]" from main.r*
export("make_interest_abundance")
make_interest_abundance <- function(data, treatment_key, dataset_name, outdir) {  #additional_taxa,
  #'*I think this is just adding the taxa for specific crops (see main.r for list)*
  # if (any(!is.na(additional_taxa))) {
  #   interest_list <- append(interest_list, additional_taxa)
  # }
  # browser()
  # interest_list <- c("Tetracyclines",
  #                    "Sulfonamides",
  #                    "Rifampin",
  #                    "Mupirocin",
  #                    "Multi-metal_resistance",
  #                    "Multi-drug_resistance",
  #                    "Multi-biocide_resistance",
  #                    "MLS",
  #                    "Iron_resistance",
  #                    "Glycopeptides",
  #                    "Drug_and_biocide_resistance",
  #                    "Drug_and_biocide_and_metal_resistance",
  #                    "Copper_resistance",
  #                    "betalactams",
  #                    "Aminoglycosides",
  #                    "Aminocoumarins")
  
  #'*original plot_data before my edits:*
  # total_ARGs = sum(data$class)
  # plot_data <- filter(data, class %in% interest_list) %>%
  #   add_other_bac() #%>%
  # #   tidy_data() %>%
  # #   calc_prop() %>%
  # #   treat_reps(treatment_key)
  browser()
  
  # plot_data_samples <- tidy_data_RA_samples(data)
  # total_ARGs = sum(plot_data_samples$Acetate_resistance + 
  #                    plot_data_samples$Acid_resistance + 
  #                    plot_data_samples$Aldehyde_resistance + 
  #                    plot_data_samples$Aminocoumarins + 
  #                    plot_data_samples$Aminoglycosides + 
  #                    plot_data_samples$Bacitracin +
  #                    plot_data_samples$betalactams + 
  #                    plot_data_samples$Biguanide_resistance + 
  #                    plot_data_samples$Biocide_and_metal_resistance +
  #                    plot_data_samples$Cationic_antimicrobial_peptides + 
  #                    plot_data_samples$Chromium_resistance + 
  #                    plot_data_samples$Copper_resistance +
  #                    plot_data_samples$Drug_and_biocide_and_metal_resistance + 
  #                    plot_data_samples$Drug_and_biocide_resistance + 
  #                    plot_data_samples$Fluoroquinolones +
  #                    plot_data_samples$Fosfomycin + 
  #                    plot_data_samples$Fusidic_acid + 
  #                    plot_data_samples$Glycopeptides +
  #                    plot_data_samples$Iron_resistance +
  #                    plot_data_samples$Mercury_resistance + 
  #                    plot_data_samples$Metronidazole + 
  #                    plot_data_samples$MLS + 
  #                    plot_data_samples$Multi_biocide_resistance + 
  #                    plot_data_samples$Multi_drug_resistance +
  #                    plot_data_samples$Multi_metal_resistance + 
  #                    plot_data_samples$Mupirocin + 
  #                    plot_data_samples$Nickel_resistance +
  #                    plot_data_samples$Nucleosides + 
  #                    plot_data_samples$Paraquat_resistance + 
  #                    plot_data_samples$Peroxide_resistance +
  #                    plot_data_samples$Phenicol + 
  #                    plot_data_samples$Phenolic_compound_resistance + 
  #                    plot_data_samples$Polyamine_resistance +
  #                    plot_data_samples$Quaternary_Ammonium_Compounds_resistance + 
  #                    plot_data_samples$Rifampin + 
  #                    plot_data_samples$Sodium_resistance +
  #                    plot_data_samples$Sulfonamides + 
  #                    plot_data_samples$Tellurium_resistance +
  #                    plot_data_samples$Tetracyclines + 
  #                    plot_data_samples$Trimethoprim + 
  #                    plot_data_samples$Zinc_resistance)
  # browser()
  plot_data_samples <- tidy_data_samples(data)  
  plot_data_samples = calc_prop_samples(plot_data_samples) %>%
    treat_reps(treatment_key)

  # plot_data_crop_av <- tidy_data_crop_av(data) %>%
  #   calc_prop_crop_av()
  #treat_reps(treatment_key)

  #' #'*sends the relative abundance data up to the ggplot function*
  plot_samples <- tidy_to_long_samples(plot_data_samples)
  # browser()
  plot_samples <- order_taxa(plot_samples)
  # browser()
  plot_samples <- plot_interest_abundance_crop_av(plot_samples)

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
  # return(plot_1)
  # utils::write.csv(plot_data_crop_av,
  #                  file = glue("{outdir}/relative_abundance_proportions_crop_av.csv"),
  #                  row.names = F)
  # browser()
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
# assumes data in tidy format and no extra coluMBs
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
  #'*as.matrix converts a data.table into a matrix, optionally using one of the coluMBs in the data.table as the matrix row names*
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
    #'*as.matrix converts a data.table into a matrix, optionally using one of the coluMBs in the data.table as the matrix row names*
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
  # browser()
  #'*at this point HBB21 sample names have the format HBB05u_t2>, while all other crops are CORe_t2_M_ON_20.  ie, have lost the year and province info*
  province_prop_data_samples <- tidy_data_province(d)
  # browser()
  province_prop_data_samples = calc_prop_province(province_prop_data_samples)

  # browser()

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
  # browser()
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
