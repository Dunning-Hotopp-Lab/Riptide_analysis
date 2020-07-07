# Evaluation of a high-throughput, cost-effective Illumina library preparation kit

Eric S. Tvedte

2020-07-07

The repository contains Supplementary Data for the manuscript, including Tables, Figures, and Files.

## Table of Contents
1. [Read representation](#read.rep)
2. [Riptide contamination analysis](#contam)
3. [Read mapping performance and genome coverage analysis](#map.cov)
4. [De novo assembly](#denovo)
5. [Metagenome analysis](#meta)
6. [System requirements](#system)
7. [Installation](#install)
8. [GitHub repository contents](#github)

### Read representation <a name="read.rep"></a>
**Raw read counts**  
Run scripts/Riptide.FileS1.Rmd using sample_data_files/raw_reads.txt

### Riptide contamination analysis <a name="contam"></a>  
**Filter 1: Species assignment, retain primary reads**
*Map reads to reference containing combined genomes from all organisms in Riptide experiment*  
bwa mem -k 23 combined.reference.fasta R1.fastq.gz R2.fastq.gz | samtools view -bho mapped.bam -  
*Sort BAM with samtools*  
samtools sort -m 2G -o sorted.bam mapped.bam  
*Mark duplicates with PicardTools*  
java -Xmx10g -jar picard.jar MarkDuplicates I=sorted.bam O=dupsmarked.bam M=dupsmetrics.txt  
*Assign reads to species*  
samtools view -F 4 -F 256 -F 1024 -F 2048 -bho primary.bam dupsmarked.bam  
samtools view primary.bam | awk -F "\t" '{print $3}' | sort -n | sed "s/\_scaf.\*//g" | uniq -c | awk '{print $2"\t"$1}' > speciesmap.txt  
Species assignment files were compiled for all libraries, with an additional row added containing numbers corresponding to the well position of each sample on the Riptide plate  
*Visualization*  
Run scripts/Riptide.FileS4 using sample_data_files/contamination.txt  

**Filter 2: Species assignment, retain primary reads with MAPQ > 20** 
samtools view -h -F 4 -F 256 -F 1024 -F 2048 -q 20 dupsmarked.bam | samtools view -bh - > filter_lowMAPQ.bam  
samtools view filter_lowMAPQ.bam | awk -F "\t" '{print $3}' | sort -n | sed "s/\_scaf.\*//g" | uniq -c | awk '{print $2"\t"$1}' > speciesmap.txt  
Compilation of data and visualization was performed as described above

### Genome coverage analysis <a name="map.cov"></a>  

**Mapped read counts**  
*Subsample read sets*  
seqkit sample -s 13 -j 4 -p #### -o R1.subsample.fastq R1.fastq.gz
-p parameter is dependent on readset  

*Map reads with bwa*  
bwa mem -k 23 reference.fasta R1.subsample.fastq R2.subsample.fastq | samtools view -bho mapped.bam  
*Sort BAM with samtools*  
samtools sort -m 2G -o sorted.bam mapped.bam  
*Mark duplicates with PicardTools*  
java -Xmx10g -jar picard.jar MarkDuplicates I=sorted.bam O=dedup.bam M=dupsmetrics.txt VALIDATION_STRINGENCY=SILENT AS=true CREATE_INDEX=true REMOVE_DUPLICATES=true
samtools view -F 4 -F 256 -F 1024 -F 2048 -bho primary.bam dedup.bam
*Visualize counts*  
Run scripts/Riptide.FileS1.Rmd using sample_data_files/mapped_reads.txt

**Histograms of genome breadth and depth of coverage**  
*Make genomebed file*  
awk -v OFS='\t' {'print $1, $2'} genome.fasta.fai > genomefile.bed  
*Make non-overlapping sliding window intervals using bedtools*  
bedtools makewindows -g genomefile.bed -w 1000 > genome_windows.bed  
*Determine genome-wide coverage using bedtools*  
bedtools coverage -a genome_windows.bed -b primary.bam -hist | grep ^all > genome_coverage.bed  
*Visualize genome-wide coverage 
Run scripts/Riptide.FileS2.Rmd - E.coli using sample_data_files/Ecoli_Riptide_depth.bed and Ecoli_KAPA_depth.bed  

**GC bias assessment**  
java -jar picard.jar CollectGcBiasMetrics I=primary.bam O=gcbias_metrics.txt CHART=gcbias_metrics.pdf S=summary_metrics.txt R=reference_sequence.fasta  
grep -v '#' gcbias_metrics.txt | awk '{print $3"\t"$4"\t"$7}' > GC.bias.final.out

*Visualize depth vs. GC*
Run scripts/Riptide.FileS3.Rmd using sample_data_files/Acinetobacter_depthvsGC.txt

### De novo assemblies <a name="denovo"></a>  
**De novo assembly of readsets with SPAdes**  
spades.py -o SPAdes_output_dir --pe1-1 R1.subsample.fastq --pe1-2 R2.subsample.fastq -m 200  

**Assembly evaluation**  
*Contiguity metrics with QUAST*  
quast.py -o $f/quast_assembly_stats -m 1000 contigs.fasta  
*Prokaryotes: conserved bacteria genes with BUSCO*  
python run_busco.py -f -c 8 -t /local/scratch/etvedte/tmp -i contigs.fasta -o busco_output_dir -l bacteria_odb9 -m geno  
*Eukaryotes: conserved metazoa genes with BUSCO*  
python run_busco.py -f -c 8 -t /local/scratch/etvedte/tmp -i contigs.fasta -o busco_output_dir -l metazoa_odb9 -m geno

### Metagenome analysis <a name="meta"></a>
**Raw read counts for metagenome libraries**
*Visualize counts*  
Run scripts/Riptide.FileS6.Rmd using sample_data_files/metagenome_data.txt  

**Assignment of LCAs to reads with kraken2**  
*Build kraken2 library*  
kraken2-build --download-taxonomy --db kraken2-nt  
kraken2-build --download library nt --db kraken2-nt  
kraken2-build --build --db kraken2-nt</span> 
*Assign reads using kraken2*  
kraken2 --threads 16 --db kraken2-nt --paired --report kraken2_report.txt  R1.fastq.gz R2.fastq.gz
*Retrieval of domain assignments for reads*  
awk '$4=="D"' kraken2_report.txt  
*Retrieval of phylum assignments for reads, minimum percentage cutoff of 0.5*  
awk '$4=="P"' kraken2_report.txt | awk '$1>0.5' -  
*Visualization*
Run scripts/Riptide.FileS6.Rmd using sample_data_files/metagenome_data.txt  

### System requirements <a name="system"></a>

R scripts were run using Windows 10 x64 with RStudio v1.1.463 using this R session:
```
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18362)
```
### Installation <a name="install"></a>
Requirements for installed packages are available in Supplementary Files. 

### GitHub repository contents <a name="github"></a>  
**supplementary_figures**  
Contains Supplementary Figures reported in the manuscript, Figure S1 - S14

**supplementary_tables**  
Contains Supplementary Tables reported in the manuscript, Table S1 - S13

**scripts**  
Contains six Rmd scripts used for data analysis:

Riptide.FileS1: Read composition, mapping performance, and indel quantification

Riptide.FileS2: Read coverage histograms

Riptide.FileS3: Contamination analysis

Riptide.FileS4: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS5: Metagenome analysis

Sample input files are provided in sample_data_files. Input and output file paths are hard-coded in the scripts, change these to run the scripts on your local system.

**htmls**  
This folder contains output html files generated from Rmd files using the R package knitr.

**sample_data_files**  
This folder contains sample input files for data analysis. 

raw+mapped_reads.txt: for use with Riptide.FileS1.Rmd

Ecoli_Riptide_depth.bed, Ecoli_KAPA_depth.bed, Shigella_KAPA_depth.bed, Shigella_KAPA_depth.bed: for use with Riptide.FileS2.Rmd

Riptide_contamination_filter1.txt: for use with Riptide.FileS3.Rmd

12 files with .coords extension: for use with Riptide.FileS4.Rmd

metagenome_data.txt: for use with Riptide.FileS5.Rmd
