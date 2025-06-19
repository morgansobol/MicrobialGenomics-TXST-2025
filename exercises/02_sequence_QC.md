# Week 2: Sequence QC

Welcome to Week 2! 
This week, you'll continue to get used to using the command line while learning how to process raw sequences and quality control (QC) them for further analysis. 

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:

- 
Preparing your raw reads from the sequencer is typically the 1st step in any genomics pipeline. Most likely you will get your raw reads back from the sequencing facility in fastq formatted files.
### What is a fastq file?

The fastq format has 4 lines per sequence: 
* the sequence identifier (header), preceded by a “@” character;
* the nucleic acid sequence itself;
* a “+” character and possibly the header information repeated;
* and the quality score information for each individual basecall, which must contain the same number of characters as letters in the sequence. 

![Breakdown of a fastq file](https://github.com/user-attachments/assets/e7a64482-ccae-4ee4-99a8-4bb83141b448)

With Illumina sequencing, the quality score information is a measure of how confident the software was when it called that particular base position. Different characters represent a specific score, i.e. here are the quality value characters in left-to-right increasing order of quality (ASCII encoding):
> !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

It's important to know that this isn’t a perfect system, as there are still confounding factors like polymerase error and other systematic errors that won’t show up in the quality score information, but performing some quality-based filtering after sequencing is essential.

We will *not* go over Demultiplexing. Demultiplexing refers to the step in processing where you’d use the barcode information in order to know which sequences came from which samples after they had all been sequenced together. Barcodes are the unique sequences that were attached to each of your invidivual samples’ genetic material before the samples got all mixed together. Demultiplexing is something that most sequence facilities will do for you nowadays. Just know that this process happens before fastq read QC. 

---

## 🧪 Exercise 1: FastQC, the most widely used tool for checking your read quality
First we need to set up the directory we are going to work in. I personally prefer to have a directory for each software we use. So let's make one for fastQC.
```bash
cd working_dir/
mkdir fastqc
cd fastqc
```
Before we run FastQC, we need to activate the Conda environment that contains the program.
```bash
conda activate seqQC
```
Now to run the program we do the following, calling the sequence reads from the rawseq folder like so:
```bash
fastqc ../../rawseqs/good_reads_R1.fq ../../fawseqs/good_reads_R2.fq
```
Now also run FastQC it on the bad_reads in the same rawseq directory. 
We should get an .html file to view the outputs.

[insert pics of good output]

[insert pics of bad output]

You can see the differences. Let's still trim both reads and see how they further improve. 

## 🧪 Exercise 2: Trimming bad reads with Trimmomatic 
Don't forget to go back one directory, into the working_dir, and create a new directory called trimmomatic.
```bash
cd ..
mkdir trimmomatic
cd trimmomatic/
```
Now we will run trimmomatic on paired-end mode, because remember from lecture, most cases reads are sequenced in the forward and reverse direction, meaning each forward read should have a paired read. 
```bash
trimmomatic PE R1.fq R2.fq xxxxxx CROP:140 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 MINLEN:140 -threads 4
```
The syntax for how to run Trimmomatic can be found in their manual (provide link), but our filtering thresholds here start with “LEADING:10”. This says cut the bases off the start of the read if their quality score is below 10, and we have the same set for the end with “TRAILING:10”. Then the sliding window parameters are 5 followed by 20, which means starting at base 1, look at a window of 5 bps and if the average quality score drops before 20, truncate the read at that position and only keep up to that point. The stringent part comes in with the MINLEN:151 at the end. Since the reads are already only 151 bps long, this means if any part of the read is truncated due to those quality metrics set above the entire read will be thrown away.

Ok, now let's check again with FastQC to see how that improved the output
```bash
cd ../fastqc/
fastqc ../trimmomatic/R1.fq R2.fq
```

