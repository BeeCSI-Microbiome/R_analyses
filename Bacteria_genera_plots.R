#Basic visualizations of kraken metagenome/pavian output


#Load in libraries
library(tidyverse)
library(cowplot)
library(dplyr)
library(ggpubr)


#read in data--------------------------------------------

#output from pavian
kraken_data <- read.csv("HBB_kraken_genera_reads.csv")

kraken_data <- tidyr::pivot_longer(
  data=kraken_data, #dataframe
  cols=5:44, # columns that we want to pivot
  names_to="Sample", # name of the column with the column titles that we are pivoting (in this case, sample names)
  values_to="reads" # name of the column with the values
)



#Add columns for environment, timepoint, replicate, and bacteria---------------------------------------------
kraken_data <- kraken_data %>%
  mutate(Environment=str_extract(Sample, 'e|u'))
kraken_data$Environment[kraken_data$Environment =="e"]<- "Exposed"
kraken_data$Environment[kraken_data$Environment =="u"]<- "Unexposed"


kraken_data <- kraken_data %>%
  mutate(Time_point=str_extract(Sample, 't1|t2|t3|t4'))



kraken_data <- kraken_data %>%
  mutate(Replicate=str_extract(Sample, '01|02|03|04|05'))


bacteria_genera <- kraken_data %>%
  mutate(bacteria=str_extract(lineage, 'bacteria'))

#Filter out bacteria
bacteria_genera <- bacteria_genera %>%
  subset(bacteria == 'bacteria')

#Transform into mean proportions/absolute abundance------------------------
# Calculate the sum of ASVs for each sample
sample_total <- 
  bacteria_genera %>%
  group_by(Sample) %>%
  mutate(sample_total = sum(reads))



# Proportion is calculated as ASV count value divided by sample total
proportions <-
  sample_total %>%
  mutate(prop = 100*reads/sample_total)


bact_gen_mean_prop <- proportions %>%
  group_by(name,Time_point, Environment) %>%
  summarise(means=mean(prop), .groups="keep",
            se= 1 * sqrt(stats::var(prop)/length(prop)
            )) %>%
  mutate(ymin = means -se, ymax = means + se)

bact_gen_mean_reads <- bacteria_genera %>%
  group_by(name,Time_point, Environment) %>%
  summarise(means=mean(reads), .groups="keep",
            se= 1 * sqrt(stats::var(reads)/length(reads)
            )) %>%
  mutate(ymin = means -se, ymax = means + se)


#functions to rearrange names and environment--------------------
#this needs to be edited to account for all bacteria
# reorder_sample_ranks <- function(level_id){
#   level_id <- factor(level_id,
#                      levels = c(
#                        "Melissococcus",
#                        "Snodgrassella",
#                        "Frischella",
#                        "Bifidobacterium",
#                        "Gilliamella",
#                        "Lactobacillus"
#                      ))
#   level_id
# }
# 
# reorder_sample_ranks_2 <- function(level_id){
#   level_id <- factor(level_id,
#                      levels = c(
#                        "Unexposed",
#                        "Exposed"
#                      ))
#   level_id
# }


#call rearrange functions------------------------------------

#bact_gen_mean_prop$name <- reorder_sample_ranks(bact_gen_mean_prop$name)
#bact_gen_mean_prop$Environment <- reorder_sample_ranks_2(bact_gen_mean_prop$Environment)

#bact_gen_mean_reads$name <- reorder_sample_ranks(bact_gen_mean_reads$name)
#bact_gen_mean_reads$Environment <- reorder_sample_ranks_2(bact_gen_mean_reads$Environment)


#exploratory analysis bacteria--------------------------------------
#stacked bacteria
ggplot(data = bact_gen_mean_prop, aes(x = Environment, y = means, fill = name)) +
  geom_bar(stat= "identity", position  = "stack") +
  #geom_jitter(aes(color = name), height = 0, width = .2) +
  labs(x = "Environment", y = "Relative abundance\n") +
  facet_grid(~ Time_point) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) 

#bacteria line chart
bar_rel<-ggplot(data = bact_gen_mean_prop, aes(x = Time_point, y = means, group = name, color = name)) +
  geom_line(linetype = "dashed", size = 1) +
  geom_point(size = 3) +
  labs(y = "Relative abundance\n")+
  facet_grid(~ Environment, scales = "free") +
  theme_bw()#+
  #theme(panel.grid.major = element_blank())

bar_rel
#bacteria line chart
line_mean_reads<- ggplot(data = bact_gen_mean_reads,  aes(x = Time_point, y = means, group = name, color = name)) +
  geom_line(linetype = "dashed", size = 1) +
  geom_point(size = 3) +
  labs(y = "Absolute abundance\n")+
  facet_grid(~ Environment, scales = "free") +
  theme_bw()#+
  #theme(panel.grid.major = element_blank())
line_mean_reads
plot_grid(bar_rel, line_mean_reads)




#boxplot for ESA
#need to subset data first
boxplot_ESA <- proportions %>%
  subset(name=='Lactobacillus' | name=='Gilliamella' | name =="Snodgrassella" | name == "Bifidobacterium" | name == "Frischella")
boxplot_ESA <- boxplot_ESA %>%
  subset(Time_point =="t1" | Time_point == "t4")
#need to reorder bacteria for plot
reorder_sample_ranks <- function(level_id){
  level_id <- factor(level_id,
                     levels = c(
                       "Lactobacillus",
                       "Gilliamella",
                       "Frischella",
                       "Snodgrassella",
                       "Bifidobacterium"
                       
                       
                       
                       
                     ))
  level_id
}
boxplot_ESA$name <- reorder_sample_ranks(boxplot_ESA$name)
boxplot_ESA$Time_point[boxplot_ESA$Time_point =="t1"]<- "T1"
boxplot_ESA$Time_point[boxplot_ESA$Time_point =="t2"]<- "T2"
boxplot_ESA$Time_point[boxplot_ESA$Time_point =="t3"]<- "T3"
boxplot_ESA$Time_point[boxplot_ESA$Time_point =="t4"]<- "T4"

ggplot(data = subset(boxplot_ESA, Environment == "Exposed"), aes(x = Time_point, y = prop, fill = name)) +
  geom_boxplot(outlier.shape = NA) +
  labs(x = "Timepoint", y = "Relative Abundance\n") +
  facet_grid(~ name, scales = "fixed") +
  #geom_jitter(aes(color = Replicate), height = 0, width = .2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        legend.position = "none")+
  scale_fill_brewer(palette = "Pastel2")+
  ylim(0,100)

