# Function scripts for various analyses and figure production including relative
# abundance, alpha diversity, and beta diversity

# Load Packages -----------------------------------------------------------
import("dplyr")
import("tidyr")
import("ggplot2")
import("vegan")
import("stats", "aggregate")
import("glue")
import("stringr")
import("openxlsx")

# Scope Variables --------------------------------------------------------
# Taxa of interest - common
interest_list <- c(
  "Lactobacillus Firm-4",
  "Lactobacillus Firm-5",
  "Other Lactobacillus",
  "Gilliamella apicola",
  "Gilliamella apis",
  "Bifidobacterium",
  "Snodgrassella alvi",
  "Frischella perrara",
  "Bartonella apis",
  "Melissococcus plutonius",
  "Paenibacillus larvae",
  "Bacteria"
)


# Wrangling Functions -----------------------------------------------------
# cleans data into tidy format
tidy_data <- function(d) {
  clean_data <- select(d, -taxRank, -taxLineage, -taxID, -depth) |>
    pivot_longer(!name, names_to = "sample", values_to = "value") |>
    pivot_wider(names_from = "name", values_from = "value")

  return(clean_data)
}

# calculates and adds Other Bacteria col
add_other_bac <- function(d) {
  column_index <- !(names(d) %in% c("taxRank", "taxLineage", "taxID", "depth", "name"))
  # Get the sum counts for taxa of interest, except Bacteria (domain), then
  # subtract their sum from Bacteria clade count. This leaves the Bacteria row
  # representing all taxa of non-interest
  df <- d[d$name != "Bacteria", column_index]
  d[d$name == "Bacteria", column_index] <-
    d[d$name == "Bacteria", column_index] - as.list(colSums(df, na.rm = T))
  # Change the grouping name
  d[d$name == "Bacteria", ]$name <- "Other Bacteria"

  return(d)
}


# calculates proportions and convert to data frame
calc_prop <- function(d) {
  sample_col <- select(d, "sample")
  prop_data <- select(d, -"sample") |>
    apply(
      MARGIN = 1,
      FUN = function(x) x / sum(x) * 100
    ) |>
    t() |>
    as.data.frame() |>
    mutate(sample_col) |>
    relocate(sample)

  return(prop_data)
}

# add in treatment and replicate cols
# Samples gathered from column names are treated with regular expressions to
# extract replicate number and treatments.
treat_reps <- function(d, treatment_key) {
  d <- d |>
    mutate(replicate = extract_replicate_string(sample),
           treatment = extract_treatment_string(sample, treatment_key))
  return(d)
}

extract_replicate_string <- function(sample_string) {
  digits <- str_extract(sample_string, "(?<=\\w{3,5})\\d\\d") |>
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
      "The following treatment code(s) detected from samples names were not found in treatment key you provided:",
      paste(missing_codes, collapse = ", ")
    ))
  }
}

# convert tidy data into long for abundance plot setup
tidy_to_long <- function(d) {
  long_data <- pivot_longer(
    data = d,
    cols = !c("sample", "treatment", "replicate"),
    names_to = "taxa",
    values_to = "value"
  )
  return(long_data)
}

# Reorders table using taxa column
order_taxa <- function(d) {
  taxa_order <- get_taxa_order(d)
  d$taxa <- factor(d$taxa, levels = taxa_order)
  d
}

# Returns a list of taxa of interest in order of descending average percent
# but with "Other Bacteria" always last
get_taxa_order <- function(d) {
  pdlong <- filter(d, taxa != "Other Bacteria")

  pd_avg <- aggregate(pdlong[, c("value")], list(pdlong$taxa), mean) |>
    arrange(value)

  taxa_order <- append(pd_avg$Group.1, "Other Bacteria")
}


# Relative Abundance ------------------------------------------------------

# plots relative abundance data for taxa of interest
# uses the taxa, value, treatment, and replicate column from data
plot_interest_abundance <- function(d) {
  color_palette <- c(
    "#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921",
    "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#508578", "#D7C1B1",
    "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0",
    "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#38333E", "#599861"
  )
  abundance_plot <- ggplot(d, aes(x = treatment, y = value, fill = taxa)) +
    scale_fill_manual(values = color_palette) +
    geom_bar(stat = "identity", colour = "black") +
    facet_grid(~replicate) +
    labs(
      title = "Relative Abundance",
      x = "Treatment",
      y = "Percentage (%)",
      fill = "Taxa"
    ) +
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    ))
}

# Returns relative abundance for taxa of interest
export("make_interest_abundance")
make_interest_abundance <- function(data, treatment_key, dataset_name,
                                    additional_taxa, outdir) {
  if (any(!is.na(additional_taxa))) {
    interest_list <- append(interest_list, additional_taxa)
  }

  plot_data <- filter(data, name %in% interest_list) |>
    add_other_bac() |>
    tidy_data() |>
    calc_prop() |>
    treat_reps(treatment_key)

  plot <- tidy_to_long(plot_data) |>
    order_taxa() |>
    plot_interest_abundance()

  ggsave(
    plot = plot,
    filename = glue("{outdir}/relative_abundance_plot.png"),
    bg = "white"
  )
  utils::write.csv(plot_data,
    file = glue("{outdir}/relative_abundance_proportions.csv"),
    row.names = F
  )
}

# separate bar plots for each treatment vs control
# CTX experiment specific, not necessary for any other experiment
export("make_separate_ctx_bars")
make_separate_ctx_bars <- function(data, treatment_key, dataset_name,
                                   additional_taxa, outdir) {
  interest_list <- append(interest_list, additional_taxa)

  clo_treat <- c("Control", "CLO")
  thi_treat <- c("Control", "THI")

  process <- filter(data, name %in% interest_list) |>
    add_other_bac() |>
    tidy_data() |>
    calc_prop() |>
    treat_reps(treatment_key)

  clo_plot <- filter(process, treatment %in% clo_treat) |>
    tidy_to_long() |>
    order_taxa() |>
    plot_interest_abundance()

  thi_plot <- filter(process, treatment %in% thi_treat) |>
    tidy_to_long() |>
    order_taxa() |>
    plot_interest_abundance()

  ggsave(plot = clo_plot, filename = glue("{outdir}/clo_abundance.png"), bg = "white")
  ggsave(plot = thi_plot, filename = glue("{outdir}/thi_abundance.png"), bg = "white")
}

# Alpha Diversity ---------------------------------------------------------
# calc alpha metrics
# assumes data in tidy format and no extra columns
calc_diversity_df <- function(d) {
  sample_col <- select(d, "sample")
  sample_data <- select(d, -"sample")
  observed_richness <- specnumber(sample_data)
  invsimpson <- diversity(sample_data, index = "invsimpson")
  simpson <- diversity(sample_data, index = "simpson")
  shannon <- diversity(sample_data, index = "shannon")
  evenness <- shannon / log(observed_richness)
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
  plot_data <- filter(d, taxRank == "S") |>
    tidy_data() |>
    calc_diversity_df() |>
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
    labs(
      title = "Alpha Diversity",
      x = "Treatment"
    )
}

# Makes all alpha div bar plots, calculates simple stats and saves them in results
export("make_all_alpha_plots")
make_all_alpha_plots <- function(data, treatment_key, dataset_name, outdir) {
  plot_data <- prep_alpha_data(data, treatment_key)

  utils::write.csv(plot_data,
    file = glue("{outdir}/alpha_div_data.csv"),
    row.names = F
  )

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
  nmds_data <- select(data, -"sample") |>
    as.matrix() |>
    metaMDS(distance = "bray") |>
    scores()
  nmds_data <- nmds_data$sites |>
    as.data.frame() |>
    mutate(sample_col) |>
    relocate(sample)

  return(nmds_data)
}

# runs anosim and saves results in a text file in results folder
calc_ano <- function(d, group_data, taxa_level, dataset_name, outdir) {
  heading <- paste(taxa_level, "ANOSIM results:\n")
  ano_data <- select(d, -"sample") |>
    as.matrix()
  rep_ano <- anosim(ano_data,
    group_data$replicate,
    distance = "bray",
    permutations = 9999
  )

  treat_ano <- anosim(ano_data,
    group_data$treatment,
    distance = "bray",
    permutations = 9999
  )

  cat(heading,
    file = glue("{outdir}/anosim.txt"), append = T
  )
  utils::capture.output(
    rep_ano,
    file = glue("{outdir}/anosim.txt"), append = T
  )
  utils::capture.output(
    treat_ano,
    file = glue("{outdir}/anosim.txt"), append = T
  )
}

# preps nmds data and runs ANOSIM on data
prep_and_ano <- function(d, treatment_key, taxa_level, dataset_name, outdir) {
  prop_data <- tidy_data(d) |>
    calc_prop()

  group_data <- prop_data |>
    calc_nmds() |>
    treat_reps(treatment_key)

  calc_ano(prop_data, group_data, taxa_level, dataset_name, outdir)

  return(group_data)
}

# plots the nmds for activity 1 data
plot_nmds_1 <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)

  hull <- data |>
    group_by(!!hull_var) |>
    slice(grDevices::chull(NMDS1, NMDS2))

  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2)) +
    geom_point(shape = 21, size = 3) +
    geom_polygon(data = hull, alpha = 0.5) +
    aes(fill = !!hull_var) +
    labs(
      title = plot_title,
      x = "NMDS1",
      y = "NMDS2",
      fill = tools::toTitleCase(h_var)
    )

  plot
}

# plots the nmds for activity 2 data
plot_nmds_2 <- function(data, h_var, plot_title) {
  hull_var <- sym(h_var)

  hull <- data |>
    group_by(!!hull_var) |>
    slice(grDevices::chull(NMDS1, NMDS2))

  plot <- ggplot(data, aes(x = NMDS1, y = NMDS2)) +
    geom_polygon(data = hull, alpha = 0.5) +
    geom_point(aes(shape = replicate), size = 3) +
    aes(fill = !!hull_var) +
    labs(
      title = plot_title,
      x = "NMDS1",
      y = "NMDS2",
      shape = "Replicate",
      fill = tools::toTitleCase(h_var)
    )

  plot
}

# plots and saves nmds' with separate hulls for replicate and treatment
act_1_nmds <- function(genus_data, speci_data, dataset_name, outdir) {
  genus_treat_nmds <- plot_nmds_1(
    genus_data,
    "treatment",
    "Genus NMDS - Brays Curtis"
  )
  genus_reps_nmds <- plot_nmds_1(
    genus_data,
    "replicate",
    "Genus NMDS - Brays Curtis"
  )
  speci_treat_nmds <- plot_nmds_1(
    speci_data,
    "treatment",
    "Species NMDS - Bray Curtis"
  )
  speci_reps_nmds <- plot_nmds_1(
    speci_data,
    "replicate",
    "Species NMDS - Bray Curtis"
  )

  ggsave(
    plot = genus_treat_nmds,
    filename = glue("{outdir}/nmds_plot_treatment_genus.png"),
    bg = "white"
  )
  ggsave(
    plot = genus_reps_nmds,
    filename = glue("{outdir}/nmds_plot_replicate_genus.png"),
    bg = "white"
  )
  ggsave(
    plot = speci_treat_nmds,
    filename = glue("{outdir}/nmds_plot_treatment_species.png"),
    bg = "white"
  )
  ggsave(
    plot = speci_reps_nmds,
    filename = glue("{outdir}/nmds_plot_replicate_species.png"),
    bg = "white"
  )
}

# plots and saves nmds' for treatment only
act_2_nmds <- function(genus_data, speci_data, dataset_name, outdir) {
  genus_treat_nmds <- plot_nmds_2(
    genus_data,
    "treatment",
    "Genus NMDS - Brays Curtis"
  )
  speci_treat_nmds <- plot_nmds_2(
    speci_data,
    "treatment",
    "Species NMDS - Bray Curtis"
  )
  ggsave(
    plot = genus_treat_nmds,
    filename = glue("{outdir}/nmds_plot_treatment_genus.png"),
    bg = "white"
  )
  ggsave(
    plot = speci_treat_nmds,
    filename = glue("{outdir}/nmds_plot_treatment_species.png"),
    bg = "white"
  )
}

# make and save nmds plots for genus and species levels using read data
# uses raw reads to calculate proportions
export("make_nmds_plots")
make_nmds_plots <- function(data, treatment_key, dataset_name, outdir) {
  file.create(glue("{outdir}/anosim.txt"))
  genus_data <- filter(data, taxRank == "G") |>
    prep_and_ano(treatment_key, "Genus", dataset_name, outdir)
  speci_data <- filter(data, taxRank == "S") |>
    prep_and_ano(treatment_key, "Species", dataset_name, outdir)

  utils::write.csv(genus_data,
    file = glue("{outdir}/nmds_plot_data_genus.csv"),
    row.names = F
  )
  utils::write.csv(speci_data,
    file = glue("{outdir}/nmds_plot_data_species.csv"),
    row.names = F
  )

  # split act 1 and 2 nmds here, based on number of treatments
  if (length(unique(genus_data$treatment)) > 2) {
    # can make convex hulls for both replicate and treatment
    act_1_nmds(genus_data, speci_data, dataset_name, outdir)
  } else {
    # can only make hulls for treatment
    act_2_nmds(genus_data, speci_data, dataset_name, outdir)
  }
}

# This code snippet may be used with the 'raw clade/taxon' table produced in
# main.R to produce a list of taxa above a certain abundance threshold.
#
# The min_abundance is the percent at which a taxon's relative abundance
# must be greater or equal to in at least n samples, where n is the value
# specified as min_samples_above_cutoff
export("taxa_cutoff_explore")
taxa_cutoff_explore <- function(taxa_count_matrix, dataset_name) {
  # Minimum percent abundance for a taxon
  min_abundance <- 1
  # Minimum number of samples that a taxon must exceed the abundance cutoff
  min_samples_above_cutoff <- 2


  pc_sc <- select(taxa_count_matrix, !matches("(taxID|taxLineage|depth)"))
  bac_sc <- pc_sc |> filter(name == "Bacteria")

  # Get proportions using Bacteria count as total
  pc_sc[c(-1, -2)] <- sapply(
    X = colnames(pc_sc[c(-1, -2)]),
    FUN = function(colname) {
      col <- pc_sc[, colname]
      col / bac_sc[[colname]] * 100
    }
  )

  taxa_of_interest_regex <- "Paenibacillus larvae|Melissococcus plutonius|Bartonella apis|Snodgrassella alvi|Lactobacillus|Frischella perrara|Bifidobacterium|Gilliamella apicola|Gilliamella apis| sp\\."

  # Get table of taxa above cutoff
  taxa_above_cutoff <- pc_sc |>
    mutate(num_samples_above_thresh = rowSums(pc_sc[c(-1, -2)] >= min_abundance)) |>
    filter(num_samples_above_thresh >= min_samples_above_cutoff &
      taxRank == "S") |>
    mutate(new_taxa_of_interest = !str_detect(name, taxa_of_interest_regex))

  # Sheet name must be substring first 31 chars of dataset_name to fit xlsx
  truncated_dataset_name <- str_sub(dataset_name, 1, 31)

  xl_file <- loadWorkbook("additional_taxa_above_1p.xlsx")
  tryCatch(
    {
      addWorksheet(xl_file, truncated_dataset_name)
    },
    error = function(cond) {
      message(cond)
    }
  )

  writeData(xl_file,
    sheet = truncated_dataset_name,
    x = select(taxa_above_cutoff, name, taxRank, new_taxa_of_interest)
  )
  saveWorkbook(xl_file, "additional_taxa_above_1p.xlsx", overwrite = TRUE)

  return(taxa_above_cutoff)
}
