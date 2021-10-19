# Microbiome analysis workflow and R analyses

This repository is intended to be the location in which we collect and share the **R scripts** for use on BeeCSI data.

Taxonomic profiling can be done with [this Snakemake pipeline](https://github.com/BeeCSI-Microbiome/taxonomic_profiling_pipeline), which produces output kraken data (taxonomic classification) and QC reports.

## A few guidelines

* A **descriptive script name** will help distinguish scripts from one another.

* **Don't worry about your R script being polished**. Even if the code not generalized, it is still useful for others for getting started. You (or others) can always update the script once it is in this repository by using _git_.

* A **header comment** at the top of the script describing what the script is for and what form the input should be (e.g. output from the Snakemake workflow, or from a different R script)

* If you'd like, you can create folders within this repository for:
    * All of your scripts (e.g. _Lance's Scripts_)
    * A series of scripts you consider to work together as a workflow
    * A single script and associated files (e.g. A README.md to describe the script and its usage, example visualizations produced by the script, etc.)

## Bioinformatics workflow

## Exploratory data analysis with Pavian

## Exploratory data analysis with Krona

## Downstream R analyses

1. Cumulative Sum Scaling normalization with `cumSum` function
2. Alpha diversity on both unnormalized and normalized data
3. Alpha 
