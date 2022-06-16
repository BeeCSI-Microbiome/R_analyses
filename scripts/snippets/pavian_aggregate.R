# This script reads in kraken reports and agglomerates them into analytic 
# matrices for both clade and taxon reads via Pavian.
#
# NO TAXA ARE FILTERED OUT AT THIS TIME
# Load packages -----------------------------------------------------------

packages <- c("dplyr",
              "purrr",
              "stringr",
              "tidyr")
lapply(packages, library, character.only = TRUE)


aggregate_reports <- function(krakenReportPaths, outpath) {
  
  krakenReportNames <-
    krakenReportPaths %>%
    map(function(x) str_replace(basename(x), "\\_report.txt$", ""))
  
  # Read in reports via Pavian
  krakenReportsPavian <-
    krakenReportPaths %>%
    map(function(x) read_report(x)) %>%
    set_names(nm = krakenReportNames)
  
  # Merge reports
  krakenReportsPavianMerged <- 
    krakenReportsPavian %>%
    map_dfr(function(x){
      x
    }, .id = "Sample")
  
  # Fix some taxa name and lineage formatting
  krakenReportsPavianMerged <- krakenReportsPavianMerged %>% 
    mutate(name = str_replace(name, "^(\\w|-)_", ""),
           taxLineage = str_replace_all(taxLineage, "\\|", ">"),
           taxLineage = str_replace_all(taxLineage, "._", ""))
  
  krakenReportsPavianMerged$taxonReads[krakenReportsPavianMerged$taxonReads == 0] <- NA
  
  merged_wide <- krakenReportsPavianMerged %>% 
    select(-percentage) %>% 
    pivot_wider(names_from = Sample, 
                values_from = c(cladeReads, taxonReads)) %>% 
    relocate(name) %>% 
    relocate(taxLineage, .after = last_col())
  
  write.csv(merged_wide, outpath, row.names = FALSE)
}

# Function taken from Pavian (https://github.com/fbreitwieser/pavian)
# Reads in kraken reports and puts them in table format
read_report <- function (myfile, has_header = NULL, check_file = FALSE) {
  first.line <- tryCatch(readLines(myfile, n = 1, warn = FALSE), 
                         error = function(e) {
                           warning("Error reading ", myfile)
                           return()
                         })
  isASCII <- function(txt) {
    if (length(txt) == 0) 
      return(FALSE)
    raw <- charToRaw(txt)
    all(raw <= as.raw(127) && (raw >= as.raw(32) | raw == 
                                 as.raw(9)))
  }
  if (length(first.line) == 0) {
    dmessage("Could not read ", myfile, ".")
    return(NULL)
  }
  tryCatch({
    if (nchar(first.line) == 0) {
      dmessage("First line of ", myfile, " is empty")
      return(NULL)
    }
  }, error = function(e) {
    dmessage(e)
    return(NULL)
  })
  if (!isTRUE(isASCII(first.line))) {
    dmessage(myfile, " is not a ASCII file")
    return(NULL)
  }
  if (is.null(has_header)) {
    has_header <- grepl("^[a-zA-Z#%\"]", first.line)
  }
  is_metaphlan3_fmt <- grepl("^#mpa_v3", first.line)
  is_metaphlan_fmt <- grepl("Metaphlan2_Analysis$", first.line)
  is_krakenu_fmt <- grepl("^.?%\treads\ttaxReads\tkmers", first.line)
  is_kaiju_fmt <- grepl("^  *%\t  *reads", first.line)
  ntabs <- lengths(regmatches(first.line, gregexpr("\t", first.line)))
  nrows <- ifelse(isTRUE(check_file), 5, -1)
  if (!is_krakenu_fmt && is_kaiju_fmt) {
    cont <- readLines(myfile)
    cont <- cont[!grepl("^-", cont)]
    cont <- sub(".*\t  *", "", cont)
    cont <- sub("; ?$", "", cont)
    report <- utils::read.delim(textConnection(cont), stringsAsFactors = FALSE)
    colnames(report) <- c("taxonReads", "taxLineage")
    report$cladeReads <- report$taxonReads
    report$taxLineage <- gsub("^", "-_", report$taxLineage)
    report$taxLineage <- gsub("; ", "|-_", report$taxLineage)
    report$taxLineage <- gsub("-_Viruses", "d_Viruses", report$taxLineage, 
                              fixed = T)
    report$taxLineage <- gsub("-_cellular organisms|-_Bacteria", 
                              "-_cellular organisms|d_Bacteria", report$taxLineage, 
                              fixed = T)
    report$taxLineage <- gsub("-_cellular organisms|-_Eukaryota", 
                              "-_cellular organisms|d_Eukaryota", report$taxLineage, 
                              fixed = T)
    report$taxLineage <- gsub("-_cellular organisms|-_Archaea", 
                              "-_cellular organisms|d_Archaea", report$taxLineage, 
                              fixed = T)
    report$taxLineage[1:(length(report$taxLineage) - 1)] <- paste0("-_root|", 
                                                                   report$taxLineage[1:(length(report$taxLineage) - 
                                                                                          1)])
    report$taxLineage[report$taxLineage == "-_unclassified"] <- "u_unclassified"
    new_counts <- integer(length = 0)
    for (j in seq_len(nrow(report))) {
      count <- report$cladeReads[j]
      tl <- report$taxLineage[j]
      tl2 <- sub("\\|[^|]*$", "", tl)
      while (tl2 != tl) {
        if (tl2 %in% names(new_counts)) {
          new_counts[tl2] <- new_counts[tl2] + count
        }
        else {
          new_counts[tl2] <- count
        }
        tl <- tl2
        tl2 <- sub("\\|[^|]*$", "", tl)
      }
    }
    report <- rbind(report, data.frame(taxonReads = 0, taxLineage = names(new_counts), 
                                       cladeReads = as.integer(new_counts)))
    tl_order <- order(report$taxLineage)
    tl_order <- c(tl_order[length(tl_order)], tl_order[-length(tl_order)])
    report <- report[tl_order, c("taxLineage", "taxonReads", 
                                 "cladeReads")]
  }
  else if (is_metaphlan3_fmt) {
    report <- tryCatch({
      utils::read.table(myfile, sep = "\t", header = F, 
                        quote = "", stringsAsFactors = FALSE, comment.char = "#", 
                        nrows = nrows, col.names = c("taxLineage", "taxID", 
                                                     "cladeReads", "additional_species"), check.names = FALSE)
    }, error = function(x) NULL, warning = function(x) NULL)
    if (is.null(report)) {
      return(NULL)
    }
    report$taxID <- sub(".*\\|", "", report$taxID)
  }
  else if (has_header) {
    report <- tryCatch({
      utils::read.table(myfile, sep = "\t", header = T, 
                        quote = "", stringsAsFactors = FALSE, comment.char = ifelse(is_metaphlan_fmt, 
                                                                                    "", "#"), nrows = nrows, check.names = FALSE)
    }, error = function(x) NULL, warning = function(x) NULL)
    if (is.null(report)) {
      return(NULL)
    }
    colnames(report)[colnames(report) %in% c("#%", "%", "clade_perc", 
                                             "perc", "percReadsClade")] <- "percentage"
    colnames(report)[colnames(report) %in% c("reads", "numReadsClade", 
                                             "n_reads_clade", "n.clade", "n-clade")] <- "cladeReads"
    colnames(report)[colnames(report) %in% c("taxReads", 
                                             "numReadsTaxon", "n_reads_taxo", "n.stay", "n-stay")] <- "taxonReads"
    colnames(report)[colnames(report) %in% c("rank", "tax_taxRank", 
                                             "level")] <- "taxRank"
    colnames(report)[colnames(report) %in% c("tax", "taxonid")] <- "taxID"
    colnames(report)[colnames(report) %in% c("indentedName", 
                                             "taxName")] <- "name"
    colnames(report)[colnames(report) %in% c("dup")] <- "kmerDuplicity"
    colnames(report)[colnames(report) %in% c("cov")] <- "kmerCoverage"
  }
  else {
    report <- NULL
    if (ntabs == 5) {
      col_names <- c("percentage", "cladeReads", "taxonReads", 
                     "taxRank", "taxID", "name")
    }
    else if (ntabs == 7) {
      col_names <- c("percentage", "cladeReads", "taxonReads", 
                     "n_unique_kmers", "n_kmers", "taxRank", "taxID", 
                     "name")
    }
    report <- tryCatch({
      utils::read.table(myfile, sep = "\t", header = F, 
                        col.names = col_names, quote = "", stringsAsFactors = FALSE, 
                        nrows = nrows)
    }, error = function(x) NULL, warning = function(x) NULL)
    if (is.null(report)) {
      dmessage(paste("Warning: File", myfile, "does not have the required format"))
      return(NULL)
    }
  }
  if (ncol(report) < 2) {
    dmessage(paste("Warning: File", myfile, "does not have the required format"))
    return(NULL)
  }
  if (is_metaphlan_fmt || is_metaphlan3_fmt) {
    colnames(report)[1] <- "taxLineage"
    colnames(report)[colnames(report) == "Metaphlan2_Analysis"] <- "cladeReads"
    report <- report[order(report$taxLineage), ]
    report$taxLineage <- gsub("_", " ", report$taxLineage)
    report$taxLineage <- gsub("  ", "_", report$taxLineage)
    report$taxLineage <- paste0("-_root|", report$taxLineage)
    root_lines <- data.frame(taxLineage = c("u_unclassified", 
                                            "-_root"), cladeReads = c(0, 100), stringsAsFactors = F)
    if ("taxID" %in% colnames(report)) {
      root_lines <- cbind(root_lines, taxID = c(0, 1))
      report <- report[, c("taxLineage", "cladeReads", 
                           "taxID")]
    }
    report <- rbind(root_lines, report)
  }
  if (all(c("name", "taxRank") %in% colnames(report)) && !"taxLineage" %in% 
      colnames(report)) {
    report$depth <- nchar(gsub("\\S.*", "", report$name))/2
    if (!all(report$depth == floor(report$depth))) {
      warning("Depth doesn't work out!")
      return(NULL)
    }
    report$name <- gsub("^ *", "", report$name)
    table(report$taxRank)
    allowed_taxRanks <- c("U", "S", "G", "F", "C", "D", "O", 
                          "K", "P")
    report$taxRank[report$taxRank == "class"] <- "C"
    report$taxRank[report$taxRank == "family"] <- "F"
    report$taxRank[report$taxRank == "genus"] <- "G"
    report$taxRank[report$taxRank == "superkingdom"] <- "D"
    report$taxRank[report$taxRank == "kingdom"] <- "K"
    report$taxRank[report$taxRank == "order"] <- "O"
    report$taxRank[report$taxRank == "phylum"] <- "P"
    report$taxRank[report$taxRank == "species"] <- "S"
    report$taxRank[report$name == "unclassified"] <- "U"
    report$taxRank[!report$taxRank %in% allowed_taxRanks] <- "-"
    report$name <- paste(tolower(report$taxRank), report$name, 
                         sep = "_")
    rownames(report) <- NULL
    report$taxLineage <- report$name
    n <- nrow(report)
    depths <- report$depth
    taxLineages <- report$name
    taxLineages_p <- as.list(seq_along(report$name))
    depth_row_tmp <- c(1:25)
    for (current_row in seq(from = 1, to = nrow(report))) {
      dcr <- depths[current_row]
      depth_row_tmp[dcr + 1] <- current_row
      if (dcr >= 1) {
        prev_pos <- depth_row_tmp[[dcr]]
        taxLineages_p[[current_row]] <- c(taxLineages_p[[prev_pos]], 
                                          current_row)
      }
    }
    report$taxLineage <- sapply(taxLineages_p, function(x) paste0(taxLineages[x], 
                                                                  collapse = "|"))
  }
  else if ("taxLineage" %in% colnames(report)) {
    taxLineages <- strsplit(report$taxLineage, "|", fixed = TRUE)
    if (!"name" %in% colnames(report)) 
      report$name <- sapply(taxLineages, function(x) x[length(x)])
    if (!"depth" %in% colnames(report)) {
      report$depth <- sapply(taxLineages, length) - 1
    }
    if (!"taxRank" %in% colnames(report)) 
      report$taxRank <- toupper(substr(report$name, 0, 
                                       1))
  }
  if (!all(c("name", "taxRank") %in% colnames(report)) || nrow(report) < 
      2 || !((report[1, "name"] == "u_unclassified" && report[2, 
                                                              "name"] == "-_root") || report[1, "name"] == "-_root")) {
    dmessage(paste("Warning: File", myfile, "does not have the required format"))
    str(report)
    return(NULL)
  }
  if (!"taxonReads" %in% colnames(report)) {
    parent <- sub("^\\(.*\\)\\|.*$", "\\1", report$taxLineage)
    taxLineages <- strsplit(report$taxLineage, "|", fixed = TRUE)
    report$taxonReads <- report$cladeReads - sapply(report$name, 
                                                    function(x) sum(report$cladeReads[parent == x]))
    report$taxonReads[report$taxonReads <= 1e-05] <- 0
  }
  report$percentage <- signif(report$cladeReads/sum(report$taxonReads), 
                              6) * 100
  if ("n_unique_kmers" %in% colnames(report)) 
    report$kmerpercentage <- round(report$n_unique_kmers/sum(report$n_unique_kmers, 
                                                             na.rm = T), 6) * 100
  if ("taxID" %in% colnames(report)) {
    std_colnames <- c("percentage", "cladeReads", "taxonReads", 
                      "taxRank", "taxID", "name")
  }
  else {
    std_colnames <- c("percentage", "cladeReads", "taxonReads", 
                      "taxRank", "name")
  }
  stopifnot(all(std_colnames %in% colnames(report)))
  report[, c(std_colnames, setdiff(colnames(report), std_colnames))]
}

krakenReportPaths <- Sys.glob("../data/ctx_2020/kraken_reports/*d0*_report.txt") %>% 
  c(Sys.glob("../data/ctx_2020/kraken_reports/*dC*_report.txt"))

outpath <- "../data/ctx_2020/CTX_dC_all_taxa.csv" 
aggregate_reports(krakenReportPaths, outpath)
