# Week 4: Isolate/Whole genome assembly & annotation

In this tutorial, we are going to assemble an unknown genome and analyze it to find out what it is and what it does.
We will be back working in the command line, not in R.


---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
- 

## 🧪 Step 1: Setting up the working environment and reading in the data
Let's make a new `conda` environment for assembly. 

```bash
conda create -n assembly -c bioconda -c conda-forge spades seqkit
```
And now let's make an environment for Anvio, follow the steps at this link:
https://anvio.org/install/macos/stable/#6-check-your-installation

Now navigate to /MicrobialGenomics-TXST-2025/data/ directory and make new directory `04_genomes`. 
```bash
cd /path/to/MicrobialGenomics-TXST-2025/data/
mkdir 04_genome
cd 04 _genome/
```

Make two directories, `data_dir` and `work_dir`. Then `cd` into `data_dir`. 
In Module 4 on Canvas, download the `raw_reads.tar.gz` where you can easily access, like downloads or your desktop.
Move this file to your current directory and unpack it. 

```bash
cp /path/to/data/raw_reads.tar.gz .
tar -xvzf raw_reads.tar.gz
ls raw_reads/
```

You should see two files: `unknown_R1.fastq.gz` and `unknown_R2.fastq.gz`. 

Now let's go back a couple of directories into `work_dir` and start processing our reads there. 
```bash
cd ../../work_dir/
mkdir fastqc spades quast
```

## 🧪 Exercise 2: Sequence QC

Like before, we are going to start with assessing our forward and reverse read quality before trimming. 
```bash
conda activate seqQC
fastqc ../../data_dir/raw_reads/*.fq -o .
```

Now run MultiQC and look at it all at once
```bash
multiqc .
open multiqc_report.html
```

We will try out a different trimming tool, _Trimmomatic_. The parameters we enter into Trimmomatic are carried out in the order in which they appear – though this is usually not the case with most programs at the command line

## 🧪 Exercise 3: Assembly


Now we’re going to put our isolate-genome assembly into the anvi’o framework and see just a few of the ways it can help us begin to look at our assembled genome.

## 🧪 Exercise 4: Anvio
For us to get our assembly into anvi’o, first we need to generate what it calls a contigs database. This contains the contigs from our assembly and information about them. The following script will organize our contigs in an anvi’o-friendly way, generate some basic stats about them, and use the program Prodigal to identify open-reading frames (ORFs).

```bash
anvi-gen-contigs-database -f contigs.fasta -o contigs.db -n unknown_genome
```
```bash
anvi-run-hmms -c contigs.db -I Bacteria_71 -T 4
```
```bash
anvi-run-hmms -c contigs.db -I Ribosomal_RNA_16S -T 4
```
```bash
anvi-setup-ncbi-cogs -T 4 # only needed the first time
anvi-run-ncbi-cogs -c contigs.db --num-threads 4
```
This will take 8-10 min.
```bash
anvi-run-kegg-kofams -c contigs.db -T 4
```
```
anvi-estimate-metabolism -c contigs.db
```

https://merenlab.org/tutorials/infant-gut/#chapter-v-metabolism-prediction

```bash
anvi-get-sequences-for-gene-calls -c contigs.db -o gene_calls.fa
```

```bash
bowtie2-build contigs.fasta contigs.btindex
bowtie2 -q -x contigs.btindex \
        -1 unknown_R1.fastq.gz \
        -2 unknown_R2.fastq.gz \
        -p 4 -S unknown_assembly.sam

samtools view -bS unknown_assembly.sam > unknown_assembly.bam

anvi-init-bam unknown_assembly.bam -o unknown_anvio.bam

anvi-profile -i unknown_anvio.bam -c contigs.db -M 1000 -T 4 --cluster-contigs -o unknown_profiled/

```

```bash
anvi-get-sequences-for-hmm-hits -c contigs.db --hmm-source Ribosomal_RNA_16S -o rRNAs.fa
```
```bash
anvi-get-sequences-for-hmm-hits -c contigs.db --hmm-sources Bacteria_71 --get-aa-sequences -o bacterial_SCGs.faa --no-wrap
```
```bash
anvi-export-functions -c contigs.db -o cog_functions.txt --annotation-sources COG20_FUNCTION,COG20_CATEGORY,COG20_PATHWAY
anvi-export-functions -c contigs.db -o kegg_functions.txt --annotation-sources KEGG_Class,KOfam,KEGG_BRITE,KEGG_Module
```

```bash
anvi-script-add-default-collection -p unknown_profiled/PROFILE.db
anvi-summarize -c contigs.db -p unknown_profiled/PROFILE.db -C DEFAULT -o unknown_assembly_summary/
```

```
anvi-interactive -c contigs.db -p unknown_profiled/PROFILE.db --title "Unknown assembly"
```

## 📝 Assignment due next class on Canvas
1. Who does your genome belong to? (Hint: what can you do with the information in rRNAs.fa file)
2. Using the provided cog_functions.txt, make a bar chart of COG functional categories. Count the number of genes per COG category letter (A, C, E, …). If a gene is annotated to multiple categories (entries separated by !!!), count it once for each category it belongs to. Label axes and add category names. Feel free to use Excel (easy) or R (advanced). If you use ChatGPT, provide a screenshot of the solution you used. 
Example from Sobol et al 2023 _BMC Genomics_
<img width="685" height="421" alt="image" src="https://github.com/user-attachments/assets/c85fa14d-da68-4267-99b9-736012eaf9ed" />









