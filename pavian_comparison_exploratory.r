# This is a script that runs using .tsv files generated from Pavian's
# comparison functionality.

# Setting working directory
# Replace "/BeeCSI/cra_data/comparisons_211006/pavian_raw" with the location
# of your .tsv files.
comparison_dir <- paste0( Sys.getenv("HOME"), "/BeeCSI/cra_data/comparisons_211006/pavian_raw")
print(comparison_dir)
getwd()
setwd(comparison_dir)
getwd()

# Loading libraries --------------------------------------------
library(tidyverse)
library(dplyr); library(plyr)

# Set some variables -------------------------------------------
num_samples = 10

# Turn all .tsv files to data frames ---------------------------
comparison_species <- read.table(file = 'comparison_species_211006.tsv', sep = '\t', header = TRUE)
comparison_phylum <- read.table(file = 'comparison_phylum_211006.tsv', sep = '\t', header = TRUE)
comparison_order <- read.table(file = 'comparison_order_211006.tsv', sep = '\t', header = TRUE)
comparison_genus <- read.table(file = 'comparison_genus_211006.tsv', sep = '\t', header = TRUE)
comparison_family <- read.table(file = 'comparison_family_211006.tsv', sep = '\t', header = TRUE)
comparison_domain <- read.table(file = 'comparison_domain_211006.tsv', sep = '\t', header = TRUE)

# Create a bar plot for phylum ---------------------------------

# Grabbing the phyla and sample name columns.
# The numbers you put here depend on what columns your file includes.
comparison_phylum_temp <- comparison_phylum[,c(1,5,6,7,8,9,10,11,12,13,14)]
# Grabbing list of sample names.
comparison_phylum_sample_names <- colnames(comparison_phylum_temp)
comparison_phylum_sample_names <- comparison_phylum_sample_names[-1]
# Shortening sample names.
comparison_phylum_sample_names <- str_remove(comparison_phylum_sample_names, ".cladeReads")
# Goal is to create a new data frame that will work with stacked bar plot.
# Samples x num_phyla
num_phyla <- nrow(comparison_phylum_temp)
samples <- c(rep(comparison_phylum_sample_names[1], num_phyla),
             rep(comparison_phylum_sample_names[2], num_phyla),
             rep(comparison_phylum_sample_names[3], num_phyla),
             rep(comparison_phylum_sample_names[4], num_phyla),
             rep(comparison_phylum_sample_names[5], num_phyla),
             rep(comparison_phylum_sample_names[6], num_phyla),
             rep(comparison_phylum_sample_names[7], num_phyla),
             rep(comparison_phylum_sample_names[8], num_phyla),
             rep(comparison_phylum_sample_names[9], num_phyla),
             rep(comparison_phylum_sample_names[10], num_phyla))
# Grabbing phyla
phyla <- comparison_phylum_temp$name
# Phyla x num_samples
phyla <- rep(c(phyla[1], phyla[2], phyla[3], phyla[4], phyla[5],
               phyla[6], phyla[7], phyla[8], phyla[9]) , num_samples)
# Populating new data frame by creating a list of all reads...
# and plopping them into the the new data frame.
phyla_cra01e <- dplyr::pull(comparison_phylum_temp, CRA01e.cladeReads)
phyla_cra01u <- dplyr::pull(comparison_phylum_temp, CRA01u.cladeReads)
phyla_cra02e <- dplyr::pull(comparison_phylum_temp, CRA02e.cladeReads)
phyla_cra02u <- dplyr::pull(comparison_phylum_temp, CRA02u.cladeReads)
phyla_cra03e <- dplyr::pull(comparison_phylum_temp, CRA03e.cladeReads)
phyla_cra03u <- dplyr::pull(comparison_phylum_temp, CRA03u.cladeReads)
phyla_cra04e <- dplyr::pull(comparison_phylum_temp, CRA04e.cladeReads)
phyla_cra04u <- dplyr::pull(comparison_phylum_temp, CRA04u.cladeReads)
phyla_cra05e <- dplyr::pull(comparison_phylum_temp, CRA05e.cladeReads)
phyla_cra05u <- dplyr::pull(comparison_phylum_temp, CRA05u.cladeReads)
phyla_values <- c(phyla_cra01e, phyla_cra01u, phyla_cra02e, 
            phyla_cra02u, phyla_cra03e, phyla_cra03u,
            phyla_cra04e, phyla_cra04u, phyla_cra05e, phyla_cra05u)
phyla_data <- data.frame(samples,phyla,phyla_values)

# Double checking to make sure my numbers aren't wack (Did I move all the numbers
# from the old data frame to the new data frame correctly)?
# If you get a false anywhere, there's a problem.
phyla_data_check <- ddply(phyla_data, .(samples), summarise, reads = sum(phyla_values))
sum(comparison_phylum$CRA01e.cladeReads) == phyla_data_check$reads[1]
sum(comparison_phylum$CRA01u.cladeReads) == phyla_data_check$reads[2]
sum(comparison_phylum$CRA02e.cladeReads) == phyla_data_check$reads[3]
sum(comparison_phylum$CRA02u.cladeReads) == phyla_data_check$reads[4]
sum(comparison_phylum$CRA03e.cladeReads) == phyla_data_check$reads[5]
sum(comparison_phylum$CRA03u.cladeReads) == phyla_data_check$reads[6]
sum(comparison_phylum$CRA04e.cladeReads) == phyla_data_check$reads[7]
sum(comparison_phylum$CRA04u.cladeReads) == phyla_data_check$reads[8]
sum(comparison_phylum$CRA05e.cladeReads) == phyla_data_check$reads[9]
sum(comparison_phylum$CRA05u.cladeReads) == phyla_data_check$reads[10]

# Plotting stacked bar plot and saving it into working directory.
# See https://www.r-graph-gallery.com/ and 
# https://ggplot2.tidyverse.org/reference/theme.html for some resources.
png('comparison_phylum.png', width = 1500,height = 1500)
ggplot(phyla_data, aes(fill=phyla, y=phyla_values, x=samples)) + 
  ggtitle("Abundance of Phyla") +
  geom_bar(position="stack", stat="identity") +
  xlab("Sample") +
  ylab("Reads") +
  theme(axis.text=element_text(size=20),
        title=element_text(size=25, face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.key.size=unit(1, "lines"),
        legend.title=element_text(size=20, face="bold"),
        legend.text=element_text(size=17))
dev.off()
