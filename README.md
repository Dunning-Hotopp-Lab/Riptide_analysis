# Evaluation of a high-throughput, cost-effective Illumina library preparation kit

Eric S. Tvedte

2019-11-08

The repository contains Supplementary Data for the manuscript, including Tables, Figures, and Files.

## System requirements

R scripts were run using Windows 10 x64 with RStudio v1.1.463 using this R session:
```
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)

## Installation
Requirements for installed packages are available in Supplementary Files. 

## (1) Scripts
This folder contains six Rmd scripts used for data analysis:

a) Riptide.FileS6: Read composition, mapping performance, and indel quantification

Riptide.FileS7: Read coverage histograms

Riptide.FileS8: Sequencing depth across GC values

Riptide.FileS9: Contamination Analysis

Riptide.FileS10: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS11: Metagenome analysis

Sample input files are provided in sample_data_files. Input and output file paths are hard-coded in the scripts, change these to run the scripts on your local system.

## (2) htmls

This folder contains output html files generated from Rmd files using the R package knitr.

## (3) sample_data_files
This folder contains sample input files for data analysis. 

Riptide.FileS6: Read composition, mapping performance, and indel quantification

Riptide.FileS7: Read coverage histograms

Riptide.FileS8: Sequencing depth across GC values

Riptide.FileS9: Contamination Analysis

Riptide.FileS10: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS11: Metagenome analysis
