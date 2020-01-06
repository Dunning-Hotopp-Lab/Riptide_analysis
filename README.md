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
6. [System requirements](#system)
7. [Installation](#install)
8. [GitHub repository contents](#github)

### Read representation analysis <a name="read.rep"></a>
**Raw read counts**  
*fastqc*  
Read pair counts doubled to determine total reads for each library.  
*Visualize counts*  
Run scripts/Riptide.FileS1.Rmd using sample_data_files/raw+mapped_reads.txt

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
Run scripts/Riptide.FileS1.Rmd using sample_data_files/raw+mapped_reads.txt

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
**Filter 1: Species assignment, retain primary reads**
*Map reads to reference containing combined genomes from all organisms in Riptide experiment*  
bwa mem -k 23 -M combined.reference.fasta R1.fastq.gz R2.fastq.gz | samtools view -bho mapped.bam  
*Sort BAM with samtools*  
samtools sort -m 2G -o sorted.bam mapped.bam  
*Mark duplicates with PicardTools*  
java -Xmx2g -jar picard.jar MarkDuplicates I=sorted.bam O=dupsmarked.bam M=dupsmetrics.txt  
*Assign reads to species*  
samtools view -F 4 -F 256 -F 1024 -F 2048 -bho primary.bam dupsmarked.bam  
samtools view primary.bam | awk -F "\t" '{print $3}' | sort -n | sed "s/\_scaf.\*//g" | uniq -c | awk '{print $2"\t"$1}' > speciesmap.txt  
Species assignment files were compiled for all libraries, with an additional row added containing numbers corresponding to the well position of each sample on the Riptide plate  
*Visualization*  
Run scripts/Riptide.FileS4 using sample_data_files/contamination.txt  

**Filter2: Species assignment, retain primary, non- multimapped reads**  
samtools view -h -F 4 -F 256 -F 1024 -F 2048 dupsmarked.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -bh - > filter_multimap.bam  
samtools view filter_multimap.bam | awk -F "\t" '{print $3}' | sort -n | sed "s/\_scaf.\*//g" | uniq -c | awk '{print $2"\t"$1}' > speciesmap.txt  
Compilation of data and visualization was performed as described above  

**Filter 3: Species assignment, retain primary reads with MAPQ > 30** 
samtools view -h -F 4 -F 256 -F 1024 -F 2048 -q 30 dupsmarked.bam | samtools view -bh - > filter_lowMAPQ.bam  
samtools view filter_lowMAPQ.bam | awk -F "\t" '{print $3}' | sort -n | sed "s/\_scaf.\*//g" | uniq -c | awk '{print $2"\t"$1}' > speciesmap.txt  
Compilation of data and visualization was performed as described above

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

**Visualizing breakpoints in E. coli assemblies**  
*Align contigs to E. coli reference using nucmer (MUMmer3)*  
nucmer --prefix NCBI Ecoli_HS_genome.fasta ecoli.spades.contigs.fasta  
delta-filter -q NCBI.delta > NCBI.filter  
show-coords -rb NCBI.filter > NCBI.filter.coords  
awk '{print $11}' NCBI.filter.coords | tail -n +6 | sort -n | uniq > NCBI.match.list  
xargs samtools faidx spades.contigs.fasta < NCBI.match.list >> matched_contigs.fasta  
*Similar alignment using only contigs with match to NCBI reference*  
nucmer --mum --minmatch 150 --prefix NCBI_match Ecoli_HS_genome.fasta matched_contigs.fasta  
delta-filter -g NCBI_match.delta > NCBI_match.filter  
show-coords -rb NCBI_match.filter > NCBI_match.coords  
cat NCBI_match.coords | tail -n +6 | awk '{print $1"\t"$2}' > NCBI_start+end.coords  
*Manual inspection of coords file to remove duplicate entries*  
awk '{print $11}' NCBI_match.coords | uniq -d  
*Visualization*  
Run scripts/Riptide.FileS.Rmd using sample_data_files/ecoli.breaks.txt

### Metagenome analysis <a name="meta"></a>
**Raw read counts for metagenome libraries**
*fastqc*  
Read pair counts doubled to determine total reads for each library.  
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

Riptide.FileS6: Read composition, mapping performance, and indel quantification

Riptide.FileS7: Read coverage histograms

Riptide.FileS8: Sequencing depth across GC values

Riptide.FileS9: Contamination Analysis

Riptide.FileS10: Visualizing breakpoints in E. coli genome assemblies

Riptide.FileS11: Metagenome analysis
