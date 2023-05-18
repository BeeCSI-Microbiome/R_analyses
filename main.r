# Main script for analysis of taxonomic classification results. This script
# encapsulates data input and calls to data processing and analysis functions,
# which are imported from accessory script files.


# 0.1 Package setup ----------------------------------------------------------
packages <- c(
  "tidyverse",
  "vegan",
  "modules",
  "data.table",
  "ggplot2",
  "glue"
)
lapply(packages, library, character.only = TRUE)
# ___________________________________________________________________________

# 0.2 Load aux scripts as modules -------------------------------------------
rsummary <- use("scripts/reads_summary.R")
ip <- use("scripts/initial_processing.R")
exploratory <- use("scripts/exploratory_functions.R")
da_ancombc <- use("scripts/da-ancombc.R")
indicsp <- use("scripts/indicator_taxa_analysis.R")
# ___________________________________________________________________________


#' Read in taxonomic count count_tables, format them, then perform analyses
#'
#' This function encapsulates taxonomic table input, processing, and subsequent
#' analyses, much of which is facilitated by calls to functions stored in
#' auxiliary scripts.
#'
#' @param dataset_name A string by which to distinguish the run of
#'   analyses/dataset. It is used exclusively in the names of output files and
#'   some figure titles.
#' @param counts_path The filepath to the aggregated count table.
#' @param treatment_key A named list which identifies the control treatment and
#'   the treatment of interest. The keys should be the sample name substrings
#'   that identify samples as belonging to the given treatment, and the values
#'   should be strings that are more descriptive. E.g. `list(u = "unexposed", e
#'   = "exposed")`
#' @param filter_string A substring identifying the treatment to be removed from
#'   the analyses (in the case of experiments with 3 treatments). The substring
#'   should uniquely identify the samples (column names in the taxa tables) with
#'   the undesired treatment. used to filter out
#' @param outdir The directory within which results will be saved.
#'
#'
#' @noRd
run_all_analyses <- function(dataset_name,
                             counts_path,
                             treatment_key,
                             filter_string = "",
                             outdir = "results") {
  # 1. Read and filter input files ------------------------------------------
  ct <- read_csv(counts_path)

  # filter the treatment that isn't being analyzed
  if (filter_string != "") {
    ct <- ct |> select(-matches(filter_string))
  }
  # _________________________________________________________________________

  # 2. Create output directories --------------------------------------------
  main_outdir <- glue("{outdir}/{dataset_name}")
  nmds_dir <- glue("{main_outdir}/nmds_anosim")
  alpha_div_dir <- glue("{main_outdir}/alpha_diversity")
  rel_abund_dir <- glue("{main_outdir}/relative_abundance")
  da_dir <- glue("{main_outdir}/differential_abundance")
  da_ancombc_dir <- glue("{da_dir}/ancombc")
  ind_sp_dir <- glue("{da_dir}/indicator_species_analysis")

  create_dir_if_nonexistant <- function(path) {
    ifelse(!dir.exists(path),
           dir.create(path, mode = "777"),
           FALSE)
  }

  lapply(
    c(outdir, main_outdir, nmds_dir, alpha_div_dir, rel_abund_dir, da_dir,
      da_ancombc_dir, ind_sp_dir),
    create_dir_if_nonexistant
  )
  # _________________________________________________________________________


  # 3. Reads summary --------------------------------------------------------
  # Perform some summary statistics before further processing
  reads_summary <- rsummary$create_summary_table(
    read_table = ct,
    dataset_name = dataset_name,
    main_outdir = main_outdir
  )
  # _________________________________________________________________________

  # 4. Formatting, filtering, and calculating clade counts ------------------
  ct <- ip$format_count_table(ct) |>
    ip$filter_table() |>
    ip$group_taxa_of_interest()

  count_tables <- ip$calculate_clade_counts(ct)

  write.csv(count_tables[["raw_clade"]],
            glue("{main_outdir}/raw_clade_counts.csv"),
            row.names = FALSE
  )
  write.csv(count_tables[["raw_taxon"]],
            glue("{main_outdir}/raw_taxon_counts.csv"),
            row.names = FALSE
  )

  # 5. ANCOMBC Differential Abundance ---------------------------------------
  # Species
  ancom_species_results <- da_ancombc$run_ancombc(
    count_table = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    rank_symbol = "S",
    outdir = da_ancombc_dir
  )
  # Genus
  ancom_genus_results <- da_ancombc$run_ancombc(
    count_table = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    rank_symbol = "G",
    outdir = da_ancombc_dir
  )

  # 6. Indicator Taxa Analysis ----------------------------------------------

  # All taxa levels combined
  indicsp$run_indicator_analysis(
    table = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    output_dir = ind_sp_dir
  )
  # Species
  indicsp$run_indicator_analysis(
    table = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    output_dir = ind_sp_dir,
    rank_symbol = "S"
  )
  # Genus
  indicsp$run_indicator_analysis(
    table = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    output_dir = ind_sp_dir,
    rank_symbol = "G"
  )

  # 7. Relative abundance ---------------------------------------------------
  # Get additional taxa above 1% rel. abundance
  taxa_cutoff_table <- exploratory$taxa_cutoff_explore(
    taxa_count_matrix = count_tables[["raw_clade"]],
    dataset_name = dataset_name
  )
  additional_taxa <- taxa_cutoff_table$name[taxa_cutoff_table$new_taxa_of_interest]

  exploratory$make_interest_abundance(
    data = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    additional_taxa = additional_taxa,
    outdir = rel_abund_dir
  )

  # 8. Alpha Diversity ------------------------------------------------------
  exploratory$make_all_alpha_plots(
    data = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    outdir = alpha_div_dir
  )

  # 9. Beta Diversity -------------------------------------------------------
  exploratory$make_nmds_plots(
    data = count_tables[["raw_clade"]],
    treatment_key = treatment_key,
    dataset_name = dataset_name,
    outdir = nmds_dir
  )

  # 10. Write the pdf report ------------------------------------------------
  rmarkdown::render(
    input = "main.qmd",
    output_file = glue("{main_outdir}/{dataset_name}_report_{Sys.Date()}.pdf"),
    params = list(main_result_dir = main_outdir,
                  dataset_name = dataset_name)
  )
}
