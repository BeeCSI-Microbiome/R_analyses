# Script for running ANCOM-BC for for pairwise differential
# abundance analysis and producing volcano plots.
# Takes two functions from exploratory_functions.R for data wrangling:
# tidy_data and treat_reps
# Takes as input: raw clade counts (should *not* be scaled/normalized)
# Uses treatment+replicate for as the model, but this can be
# customized in the ancombc function itself.

# Author(s): Jonathan Ho
# Updated: Apr 14, 2022

modules::import("tidyr")
modules::import("stringr")
modules::import("ANCOMBC")
modules::import("microbiome")
modules::import("dplyr")
modules::import("tibble")
modules::import("ggplot2")
modules::import("utils")
modules::import("glue")
modules::import("phyloseq")


taxa_lvl_key <- c(D="Domain",
                  P="Phylum",
                  C="Class",
                  O="Order",
                  F="Family",
                  G="Genera",
                  S="Species")

# User Defined Variables --------------------------------------------------
base_title <- 'Cage Control vs'


# Functions from exploratory_functions.R ----------------------------------

# cleans data into tidy format
tidy_data <- function(d) {
  clean_data <- select(d,-taxRank,-taxLineage, -taxID, -depth) %>%
    pivot_longer(!name, names_to = "sample", values_to = "value") %>%
    pivot_wider(names_from = "name", values_from = "value")

  return(clean_data)
}

# add in treatment and replicate cols
# Samples gathered from column names are treated with regular expressions to
# extract replicate number and treatments.
treat_reps <- function(d, treatment_key) {
  d <- d %>% mutate(replicate = extract_replicate_string(sample),
                    treatment = extract_treatment_string(sample, treatment_key))
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
      "The following treatment code(s) detected from samples names were not found in treatment key you provided:",
      paste(missing_codes, collapse = ", ")))
  }
}


# Functions for ANCOM-BC --------------------------------------------------
# setup data and object for ancombc
prep_data <- function(count_table, treatment_key, rank, meta_table) {
  data <- count_table %>%
    filter(taxRank == rank)

  abun_data <- select(data,-taxRank,-taxLineage, -taxID, -depth) %>%
    column_to_rownames('name') %>%
    as.matrix()

  metadata <- tidy_data(data) %>%
    select(sample) %>%
    treat_reps(treatment_key)

  if (!is.null(meta_table)) {
    metadata <- metadata %>%
      merge(select(meta_table, Sample, Lab), by.x = "sample", by.y = "Sample")
  }
  metadata <- metadata %>% column_to_rownames('sample')

  OTU <- otu_table(abun_data, taxa_are_rows = T)
  samples <- sample_data(metadata)

  phylo_obj <- phyloseq(OTU, samples)

  sample_data(phylo_obj)$treatment <-
    as.factor(sample_data(phylo_obj)$treatment) %>%
    stats::relevel(treatment_key[[1]])

  return(phylo_obj)
}

# visualizes and saves plots and csv
visualize_save <-
  function(mod,
           treatment_key,
           dataset_name,
           rank,
           outdir) {
    # only get treatment related log-fold-changes (lfc), adj p-vals, and diffs
    df <- data.frame(mod$res) %>%
      column_to_rownames(var = "taxon") %>%
      select(contains('treatment')) %>%
      select(starts_with('lfc_') | starts_with('p_') | starts_with('q_') |
               starts_with('diff_')) %>%
      mutate(across(
        starts_with('lfc_'),
        .fns = function(x)
          log2(exp(x)),
        .names = "log2_{col}"
      )) %>%
      relocate(starts_with("log2_"))

    # loop through treatment names and do everything inside
    controls <- c('control', 'unexposed')
    all_treatments <- unlist(treatment_key, use.names = FALSE)
    just_treats <- all_treatments[!all_treatments %in% controls]

    for (i in just_treats) {
      # setup identifiable strings for i-th treatment
      diff_abun_i <- paste('diff_abun_', i, sep = '') %>%
        sym()
      q_val_i <- paste('q_treatment', i, sep = '') %>%
        sym()
      log2_lfc_i <- paste('log2_lfc_treatment', i, sep = '') %>%
        sym()
      plot_title_i <-
        paste(base_title, i, '-', str_to_title(rank), '-', dataset_name)

      # determines whether DA are increases or decreases
      df <- mutate(df,!!diff_abun_i := ifelse(
        df[, grep(q_val_i,
                  colnames(df))] < 0.05,
        ifelse(df[, grep(log2_lfc_i,
                         colnames(df))] >= 0,
               'Increase',
               'Decrease'),
        'No Change'
      ))

      # adjust factor levels for plotting
      df[, grep(diff_abun_i, colnames(df))] <-
        factor(df[, grep(diff_abun_i, colnames(df))],
               levels = c('Increase', 'No Change', 'Decrease'))

      # plot and save
      da_plot <-
        volc_plot(df, log2_lfc_i, q_val_i, diff_abun_i, plot_title_i)
      ggsave(da_plot,
             filename = glue("{outdir}/{plot_title_i}.png"))
    }
    # Write the full ancom results
    write.csv(
      arrange(df, row.names(df)),
      file = glue("{outdir}/{dataset_name}_{rank}_ancombc.csv"),
      row.names = T
    )
    # Format and write the collaborator-formatted ancom results
    df <- df %>%
      rownames_to_column(var = "taxa_name") %>%
      select(matches("^(taxa_name|p_|q_|log2).*")) %>%
      rename_with(.fn = ~ "log2_fold_change", .cols = starts_with("log2")) %>%
      rename_with(.fn = ~ "p-value", .cols = starts_with("p_")) %>%
      rename_with(.fn = ~ "corrected_p-value", .cols = starts_with("q_"))

    write.csv(
      df,
      file = glue("{outdir}/{dataset_name}_{rank}_ancombc_collaborator.csv"),
      row.names = FALSE
    )
  }

# volcano plot helper
volc_plot <-
  function(df,
           log2_lfc_i,
           q_val_i,
           diff_abun_i,
           plot_title_i) {
    # determine x range
    x_default = 1.35
    x_data = max(abs(df[, grep(log2_lfc_i, colnames(df))]))
    x_use = max(x_default, x_data) + 0.15

    # determine y range
    y_default = 1.40
    y_data = -log10(min(df[, grep(q_val_i, colnames(df))]))
    y_use = max(y_default, y_data) + 0.10

    # plot
    da_plot <- ggplot(df,
                      aes(
                        x = !!log2_lfc_i,
                        y = -log10(!!q_val_i),
                        colour = !!diff_abun_i
                      )) +
      geom_point(alpha = 0.5, size = 3) +
      scale_color_manual(values = c(
        'Increase' = 'blue',
        'No Change' = 'black',
        'Decrease' = 'red'
      )) +
      xlim(0 - x_use,
           0 + x_use) +
      ylim(0, y_use) +
      geom_hline(
        yintercept = -log10(0.05),
        lty = 4,
        col = "black",
        lwd = 0.8
      ) +
      labs(x = "log2(fold change)",
           y = "-log10(adj. p-value)",
           title = plot_title_i) +
      theme(legend.position = "right",
            legend.title = element_blank())

    return(da_plot)
  }


export("run_ancombc")
run_ancombc <- function(count_table,
                        treatment_key,
                        dataset_name,
                        rank_symbol,
                        outdir,
                        meta_table) {

  has_multiple_treatments <- length(treatment_key) > 2
  phylo_obj <-
    prep_data(count_table, treatment_key, rank_symbol, meta_table)

  if (has_multiple_treatments) {
    ancom_result <- ancombc2(data = phylo_obj,
                             fix_formula = 'treatment',
                             rand_formula = '(1 | replicate)',
                             p_adj_method = "BH",
                             group = "treatment",
                             mdfdr_control = list(fwer_ctrl_method = "BH", B = 100),
                             dunnet = TRUE,
                             global = TRUE)
  } else {
    if (length(unique(sample_data(phylo_obj)$Lab)) > 1) {
      rand_formula_str <- '(1 | replicate) + (1 | Lab)'
    } else {
      rand_formula_str <- '(1 | replicate)'
    }
    ancom_result <- ancombc2(data = phylo_obj,
                             fix_formula = 'treatment',
                             rand_formula = rand_formula_str,
                             p_adj_method = "BH")
  }

  visualize_save(ancom_result, treatment_key, dataset_name, taxa_lvl_key[rank_symbol], outdir)
  return(ancom_result)
}
