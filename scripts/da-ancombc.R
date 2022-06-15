# Script for running ANCOM-BC for for pairwise differential
# abundance analysis and producing volcano plots.
# Takes two functions from exploratory_functions.R for data wrangling:
# tidy_data and treat_reps
# Takes as input: raw clade counts (should *not* be scaled/normalized)
# Uses treatment+replicate for as the model, but this can be
# customized in the ancombc function itself.

# Author(s): Jonathan Ho
# Updated: Apr 14, 2022

import("tidyr")
import("stringr")
import("ANCOMBC")
import("microbiome")
import("dplyr")
import("tibble")
import("ggplot2")
import("utils")
import("glue")
import("phyloseq")

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
  } else if (all(str_detect(sample_string, "(?<=[[:upper]]{3}\\d\\d)(e|u)"))) {
    # Or activity 2 pattern?
    treatment_codes <-
      str_extract(sample_string, "(?<=[[:upper]]{3}\\d\\d)(e|u)")
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
prep_data <- function(count_table, treatment_key, rank) {
  data <- count_table %>%
    filter(taxRank == rank)
  
  abun_data <- select(data,-taxRank,-taxLineage, -taxID, -depth) %>%
    column_to_rownames('name') %>%
    as.matrix()
  
  metadata <- tidy_data(data) %>%
    select(sample) %>%
    treat_reps(treatment_key) %>%
    column_to_rownames('sample')

  OTU <- otu_table(abun_data, taxa_are_rows = T)
  samples <- sample_data(metadata)
  
  phylo_obj <- phyloseq(OTU, samples)
  return(phylo_obj)
}

# visualizes and saves plots and csv
visualize_save <-
  function(mod,
           treatment_key,
           dataset_name,
           rank,
           outdir) {
    # only get treatment related betas, adj p-vals, and diffs
    df <- data.frame(mod$res) %>%
      select(contains('treatment')) %>%
      select(contains('beta') | contains('q_val') |
               contains('diff')) %>%
      mutate(across(
        contains('beta'),
        .fns = function(x)
          log2(exp(x)),
        .names = "log2_{col}"
      )) %>%
      relocate(contains("log2"))
    
    # loop through treatment names and do everything inside
    controls <- c('control', 'unexposed')
    all_treatments <- unlist(treatment_key, use.names = FALSE)
    just_treats <- all_treatments[!all_treatments %in% controls]

    for (i in just_treats) {
      # setup identifiable strings for i-th treatment
      diff_abun_i <- paste('diff_abun_', i, sep = '') %>%
        sym()
      q_val_i <- paste('q_val.treatment', i, sep = '') %>%
        sym()
      beta_i <- paste('log2_beta.treatment', i, sep = '') %>%
        sym()
      plot_title_i <-
        paste(base_title, i, '-', str_to_title(rank), '-', dataset_name)
      
      # determines whether DA are increases or decreases
      df <- mutate(df,!!diff_abun_i := ifelse(
        df[, grep(q_val_i,
                  colnames(df))] < 0.05,
        ifelse(df[, grep(beta_i,
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
        volc_plot(df, beta_i, q_val_i, diff_abun_i, plot_title_i)

      ggsave(da_plot,
             filename = glue("{outdir}/{plot_title_i}.png"))
    }
    write.csv(
      df,
      file = glue("{outdir}/{dataset_name}_{rank}_ancombc.csv"),
      row.names = T
    )
    
  }

# volcano plot helper
volc_plot <-
  function(df,
           beta_i,
           q_val_i,
           diff_abun_i,
           plot_title_i) {
    # determine x range
    x_default = 1.35
    x_data = max(abs(df[, grep(beta_i, colnames(df))]))
    x_use = max(x_default, x_data) + 0.15
    
    # determine y range
    y_default = 1.40
    y_data = -log10(min(df[, grep(q_val_i, colnames(df))]))
    y_use = max(y_default, y_data) + 0.10
    
    # plot
    da_plot <- ggplot(df,
                      aes(
                        x = !!beta_i,
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
      geom_vline(
        xintercept = c(-1, 1),
        lty = 4,
        col = "black",
        lwd = 0.8
      ) +
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
                        outdir) {
  phylo_obj <-
    prep_data(count_table, treatment_key, rank_symbol)
  ancom_result <- ancombc(phylo_obj,
                          formula = 'treatment+replicate',
                          p_adj_method = "BH")
  visualize_save(ancom_result, treatment_key, dataset_name, rank_symbol, outdir)
}
