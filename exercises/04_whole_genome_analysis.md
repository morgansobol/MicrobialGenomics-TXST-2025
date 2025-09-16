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
conda install -c conda-forge mamba --yes
mamba create -n assembly -c bioconda -c conda-forge spades quast 
```

Now navigate to /MicrobialGenomics-TXST-2025/data/ directory and make new directory `04_genomes`. 
```bash
cd /path/to/MicrobialGenomics-TXST-2025/data/
mkdir 04_genomes
cd 04 _genomes/
```

Make two directories, `data_dir` and `work_dir`. Then `cd` into `data_dir`. 
In Module 4 on Canvas, download the `raw_reads.tar.gz` where you can easily access, like downloads or your desktop.
Move this file to your current directory and unpack it. 

```bash
cp /path/to/data/raw_reads.tar.gz .
tar -xvzf raw_reads.tar.gz
ls raw_reads/
```

You should see two files: `unknwon_R1.fastq.gz` and `unknwon_R2.fastq.gz`. 

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



