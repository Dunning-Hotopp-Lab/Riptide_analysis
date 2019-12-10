# Evaluation of a high-throughput, cost-effective Illumina library preparation kit

Eric S. Tvedte

2019-11-08

The repository contains Supplementary Data for the manuscript, including Tables, Figures, and Files.

## Table of Contents
1. [Read representation analysis](#read.rep)
2. [Genome coverage analysis](#cov)
3. [Riptide contamination analysis](#contam)
4. [De novo assemblies](#denovo)
5. [Metagenome analysis](#meta)

### Read representation analysis <a name="read.rep"></a>
**Raw read counts**  
*fastqc*  
Read pair counts doubled to determine total reads for each library.  
*Visualize counts*  
Run scripts/FileS1 - raw read counts using sample_data_files/raw+mapped_reads.txt

**Mapped read counts**  
*Map reads with bwa*  
bwa mem -k 23 -M reference.fasta R1.fastq.gz R2.fastq.gz | samtools view -bho mapped.bam  
*Sort BAM with samtools*  
samtools sort -m 2G -o sorted.bam mapped.bam  
*Mark duplicates with PicardTools*  
java -Xmx2g -jar picard.jar MarkDuplicates I=sorted.bam O=dupsmarked.bam M=dupsmetrics.txt  
*Count reads, filtering unmapped and duplicate*  
samtools view -F 4 -F 1024 -c dupsmarked.bam  
*Visualize counts*  
Run scripts/Riptide.FileS1.Rmd - mapped read counts using sample_data_files/raw+mapped_reads.txt 
*Count reads, filtering unmapped, duplicate, secondary*  
samtools view -F 4 -F 256 -F 1024 -F 2048 -c dupsmarked.bam  
*Visualize counts*  
Run scripts/Riptide.FileS1.Rmd - primary read counts using sample_data_files/raw+mapped_reads.txt

**Quantification of indels in read libraries**
*Generate filtered BAM to retain primary reads*  
samtools view -F 4 -F 256 -F 1024 -F 2048 -bho primary.bam dupsmarked.bam   
*Generate pileup of mapped reads using samtools*  
samtools mpileup -R -d 0 -f reference.fasta primaryreads.bam > mpileup.txt
*Calculate insertions in pileup*  
cat mpileup.txt | awk '{print $5}' | grep -o '\+[0-9]\+[ACGTNacgtn]\+' | wc -l > insertion_counts.txt
*Calculate deletions in pileup* 
cat mpileup.txt | awk '{print $5}' | grep -o '\-[0-9]\+[ACGTNacgtn]\+' | wc -l > deletion_counts.txt
*Truncate Theileria 2x300 NEBNext library*  
fastx_trimmer -l 150 -i Tp_NEBNExt_fwd.fastq -o Tp_NEBNExt_trimmed_fwd.fastq  
fastx_trimmer -l 150 -i Tp_NEBNExt_rev.fastq -o Tp_NEBNExt_trimmed_rev.fastq  
Perform mapping, pileup, calculations as described above.

### Genome coverage analysis <a name="cov"></a>  
**Histograms of genome breadth and depth of coverage**  
*Generate filtered BAM to retain primary reads*  
samtools view -F 4 -F 256 -F 1024 -F 2048 -bho primary.bam dupsmarked.bam  
*Merge readsets in Aspergillus, Brugia, Plasmodium*  
samtools merge merge.Aspergillus.bam Aspergillus.1.bam Aspergillus.2.bam Aspergillus.3.bam ...  
*Make genomebed file*  
awk -v OFS='\t' {'print $1, $2'} genome.fasta.fai > genomefile.bed  
*Make non-overlapping sliding window intervals using bedtools*  
bedtools makewindows -g genomefile.bed -w 1000 > genome_windows.bed  
*Determine genome-wide coverage using bedtools*  
bedtools coverage -a genome_windows.bed -b primary.bam -hist | grep ^all > genome_coverage.bed  
*Visualize genome-wide coverage: continuous y-axis*  
Run scripts/Riptide.FileS2.Rmd - Klebsiella using sample_data_files/Klebsiella_coverage.txt  
*Visualize genome-wide coverage: discontinuous y-axis*  
Run scripts/Riptide.FileS2.Rmd - Plasmodium using sample_data_files/Plasmodium_coverage.txt  
*Make overlapping sliding window intervals using bedtools*  
bedtools makewindows -g genomefile.bed -w 1000 -s 500 > genome_windows_ovl.bed  

**Sequencing depth versus genome GC content**  
*Determine GC content of sliding windows*  
bedtools nuc -fi genome.fasta -bed genome_windows_ovl.bed | awk '{print $1,$2, $3, $5}' - | tail -n +2 - > GC.txt  
*Determine mode depth for each sliding window*  
bedtools coverage -sorted -d -a genome_windows_ovl.bed -b primary.bam| sort -k1,1 -k2,2n | groupBy -g 1,2,3 -c 5 -o mode | awk '{print $4}' - > mode.txt  
*Organize output files*  
paste GC.txt mode.txt > data.txt  
echo -e "scaffold\twindow_start\twindow_end\tGC_content\tmode_depth" > header.txt  
cat header.txt data.txt | sed -e "s/ /\t/g" > depth+GCstats.txt  
*Visualize depth vs. GC*
Run scripts/Riptide.FileS3.Rmd using sample_data_files/Acinetobacter_depthvsGC.txt

### Riptide contamination analysis <a name="contam"></a>


### De novo assemblies <a name="denovo"></a>
**De novo assembly of readsets with SPAdes**  
*Single paired end datasets*  
spades.py -o SPAdes_output_dir --pe1-1 fwd.fastq --pe1-2 rev.fastq -m 200  
*Merged datasets with multiple paired end files*  
spades.py -o SPAdes_output_dir --dataset merged.readsets.yaml -m 200  

**Assembly evaluation**  
*Contiguity metrics with QUAST*  
quast.py -o $f/quast_assembly_stats -m 1000 contigs.fasta  
*Prokaryotes: conserved bacteria genes with BUSCO*  
python run_busco.py -f -c 8 -t /local/scratch/etvedte/tmp -i contigs.fasta -o busco_output_dir -l bacteria_odb9 -m geno  
*Eukaryotes: conserved metazoa genes with BUSCO*  
python run_busco.py -f -c 8 -t /local/scratch/etvedte/tmp -i contigs.fasta -o busco_output_dir -l metazoa_odb9 -m geno

**Subsampling non-Riptide libraries**
*Retrieve 75% random sample of readsets using seqtk*  
seqtk sample -s 13 fwd.fastq 0.75 > fwd.0.75.fastq  
seqtk sample -s 13 rev.fastq 0.75 > rev.0.75.fastq  
Repeat to retrieve 50% and 25% of non-Riptide library reads.  
Perform spades and QUAST as above with subsampled datasets.

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
This folder contains Supplementary Tables reported in the manuscript, Table S1 - SX

## 3. scripts
This folder contains six Rmd scripts used for data analysis:

Riptide.FileS6: Read composition, mapping performance, and indel quantification

Riptide.FileS7: Read coverage histograms

Riptide.FileS8: Sequencing depth across GC values

Riptide.FileS9: Contamination Analysis

Riptide.FileS10: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS11: Metagenome analysis

Sample input files are provided in sample_data_files. Input and output file paths are hard-coded in the scripts, change these to run the scripts on your local system.

## 4. htmls

This folder contains output html files generated from Rmd files using the R package knitr.

## 5. sample_data_files
This folder contains sample input files for data analysis. 

Riptide.FileS6: Read composition, mapping performance, and indel quantification

Riptide.FileS7: Read coverage histograms

Riptide.FileS8: Sequencing depth across GC values

Riptide.FileS9: Contamination Analysis

Riptide.FileS10: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS11: Metagenome analysis
