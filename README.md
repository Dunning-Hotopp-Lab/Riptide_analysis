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
```
## Installation
Requirements for installed packages are available in Supplementary Files. 

## 1. supplementary_figures
This folder contains Supplementary Figures reported in the manuscript, Figure S1 - S12

## 2. supplementary_tables
This folder contains Supplementary Tables reported in the manuscript, Table S1 - S3

## 3. supplementary_files
This folder contains Supplementary Files reported in the manuscript, File S1 - S5

## 4. scripts
This folder contains six Rmd scripts used for data analysis:

Riptide.FileS6: Read composition, mapping performance, and indel quantification

Riptide.FileS7: Read coverage histograms

Riptide.FileS8: Sequencing depth across GC values

Riptide.FileS9: Contamination Analysis

Riptide.FileS10: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS11: Metagenome analysis

Sample input files are provided in sample_data_files. Input and output file paths are hard-coded in the scripts, change these to run the scripts on your local system.

## 5. htmls

This folder contains output html files generated from Rmd files using the R package knitr.

## 6. sample_data_files
This folder contains sample input files for data analysis. 

Riptide.FileS6: Read composition, mapping performance, and indel quantification

Riptide.FileS7: Read coverage histograms

Riptide.FileS8: Sequencing depth across GC values

Riptide.FileS9: Contamination Analysis

Riptide.FileS10: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS11: Metagenome analysis
