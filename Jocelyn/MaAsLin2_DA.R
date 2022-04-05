if(!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("Maaslin2")

#install.packages("devtools")


#Load MaAsLin 2 package into the R environment
import("devtools")
# install_github("biobakery/Maaslin2")
import("Maaslin2")
import("ggplot2")
import("glue")
import("utils")
# Maaslin2

export("MaAsLin2_DA_function")
MaAsLin2_DA_function <- function(kraken_matrix_dir, mas_metadata_filepath, mas_dir){


#{kraken_matrix_dir}
clade_input <- Sys.glob(glue("{kraken_matrix_dir}/krakenAnalytical_cladeReads.tsv"))
taxon_input <- Sys.glob(glue("{kraken_matrix_dir}/krakenAnalytical_taxonReads.tsv"))

clade_input_table <- read.table(clade_input, sep="\t", header=T, row.names = 1, 
                                comment.char = "", quote="", check.names = F)
taxon_input_table <- read.table(taxon_input, sep="\t", header=T, row.names = 1, 
                                comment.char = "", quote="", check.names = F)
#metadata files-----
metadata_txt = Sys.glob(mas_metadata_filepath)

clade_groupings <- read.table(metadata_txt, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
taxon_groupings <- read.table(metadata_txt, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)

#output files-----
output_file = Sys.glob(mas_dir)


#number of samples--------------------------------------------------------------------
#clade-----
clade_sample_num <- length(colnames(clade_input_table))
clade_grouping_num <- length(rownames(clade_groupings))

if(clade_sample_num != clade_grouping_num){
  message("The number of samples in the clade_input_table and the metadata table are unequal")
  message("Will remove any samples that are not found in either the clade_input_table or the metadata table")
}

#taxon-----
taxon_sample_num <- length(colnames(taxon_input_table))
taxon_grouping_num <- length(rownames(taxon_groupings))

if(taxon_sample_num != taxon_grouping_num){
  message("The number of samples in the taxon_input_table and the metadata table are unequal")
  message("Will remove any samples that are not found in either the taxon_input_table or the metadata table")
}

#reformatting input tables-------------------------------------------------------------------------------------
#This is rearranging the data and metadata to a format with the samples in the same order.  But why?  The MaAsLin2 README says that's not necessary.
#   Maybe because they converted it to a data.table? That messes with the Maaslin2 function?

#clade-----
if(identical(colnames(clade_input_table), rownames(clade_groupings))==T){
  message("Metadata and clade_input_table are in the same order")
}else{
  rows_to_keep <- intersect(colnames(clade_input_table), rownames(clade_groupings))
  clade_groupings <- clade_groupings[rows_to_keep,,drop=F]
  clade_input_table <- clade_input_table[,rows_to_keep]
  if(identical(colnames(clade_input_table), rownames(clade_groupings))==T){
    message("Metadata table was re-arrange to be in the same order as the clade_input_table")
    message("A total of ", clade_sample_num-length(colnames(clade_input_table)), " from the clade_input_table")
    message("A total of ", clade_grouping_num-length(rownames(clade_groupings)), " from the metadata table")
  }else{
    stop("Unable to match samples between the clade_input_table and metadata table")
  }
}

#taxon---------
if(identical(colnames(taxon_input_table), rownames(taxon_groupings))==T){
  message("Metadata and taxon_input_table are in the same order")
}else{
  rows_to_keep <- intersect(colnames(taxon_input_table), rownames(taxon_groupings))
  taxon_groupings <- taxon_groupings[rows_to_keep,,drop=F]
  taxon_input_table <- taxon_input_table[,rows_to_keep]
  if(identical(colnames(taxon_input_table), rownames(taxon_groupings))==T){
    message("Metadata table was re-arrange to be in the same order as the taxon_input_table")
    message("A total of ", taxon_sample_num-length(colnames(taxon_input_table)), " from the taxon_input_table")
    message("A total of ", taxon_grouping_num-length(rownames(taxon_groupings)), " from the metadata table")
  }else{
    stop("Unable to match samples between the taxon_input_table and metadata table")
  }
}


library(Maaslin2)

clade_input_table <- data.frame(t(clade_input_table), check.rows = F, check.names = F, stringsAsFactors = F)
taxon_input_table <- data.frame(t(taxon_input_table), check.rows = F, check.names = F, stringsAsFactors = F)

#call the function (Nearing's call format)-------------------------------------------------------------------------------
#clade-----
fit_data <- Maaslin2(
  clade_input_table, clade_groupings, paste(output_file, "/", "clade_results", sep = "", collapse = ""), transform = "AST",
  fixed_effects = c(colnames(clade_groupings[1])),
  random_effects = c(colnames(clade_groupings[2])),
  standardize = FALSE, plot_heatmap = F, plot_scatter = F)

#taxon----
fit_data <- Maaslin2(
  taxon_input_table, taxon_groupings, paste(output_file, "/", "taxon_results", sep = "", collapse = ""), transform = "AST",
  fixed_effects = c(colnames(taxon_groupings[1])),
  random_effects = c(colnames(taxon_groupings[2])),
  standardize = FALSE, plot_heatmap = F, plot_scatter = F)

#Convert results to .csv-------------------------------------------------------------------------------------------
#clade---
#significant results
table_sig_clade = read.table(paste(output_file, "/", "clade_results", "/", "significant_results.tsv", sep = "", collapse = ""))
write.csv(table_sig_clade,
            file = paste(output_file, "/", "clade_results", "/", "clade_significant_results.csv", sep = "", collapse = ""),
            quote = F, row.names = F)
#all results
table_all_clade = read.table(paste(output_file, "/", "clade_results", "/", "all_results.tsv", sep = "", collapse = ""))
write.csv(table_all_clade,
          file = paste(output_file, "/", "clade_results", "/", "clade_all_results.csv", sep = "", collapse = ""),
          quote = F, row.names = F)

#taxon---
#significant results
table_sig_taxon = read.table(paste(output_file, "/", "taxon_results", "/", "significant_results.tsv", sep = "", collapse = ""))
write.csv(table_sig_taxon,
          file = paste(output_file, "/", "taxon_results", "/", "taxon_significant_results.csv", sep = "", collapse = ""),
          quote = F, row.names = F)
#all results
table_all_taxon = read.table(paste(output_file, "/", "taxon_results", "/", "all_results.tsv", sep = "", collapse = ""))
write.csv(table_all_taxon,
          file = paste(output_file, "/", "taxon_results", "/", "taxon_all_results.csv", sep = "", collapse = ""),
          quote = F, row.names = F)



sessionInfo()


}
