# Week 2: Sequence QC

Welcome to Week 2! 
This week, you'll continue to get used to using the command line while learning how to process raw sequences and quality control (QC) them for further analysis. 

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:

- Understand the fastq sequence format.
- Evaluate raw sequencing data quality and explain how filtering/trimming decisions affect read retention.
- Learn the importance of read QC for downstream analyses. 
  
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

## 🧪 Exercise 1: FastQC, checking your read quality
First we need to set up the directory we are going to work in. Head into the MicrobialGenomics-TXST-2025/data/02_sequenceQC/ directory in your MicrobialGenomics-TXST-2025 from last week.
```bash
 cd MicrobialGenomics-TXST-2025/data/02_sequenceQC/
```

Now let's download the data we will work with
```bash
wget https://raw.githubusercontent.com/morgansobol/MicrobialGenomics-TXST-2025/main/data/02_sequenceQC/data_dir.tar.gz
```
If `wget` does not work, try curl instead:
```bash
curl -L -O https://raw.githubusercontent.com/morgansobol/MicrobialGenomics-TXST-2025/main/data/02_sequenceQC/data_dir.tar.gz
```

Now this is a compressed file, also called a "tar ball". To unpack it, like opening a zip file, we run:
```bash
tar -xzvf data_dir.tar.gz
```

This should have unpacked a new directory called data_dir.

Let's set up the rest of our environment to process the data.

I personally prefer to have a directory for each software we use. So let's make one for fastQC.
```bash
mkdir working_dir/
cd working_dir/
mkdir fastqc
cd fastqc
```

Before we run FastQC, we need to install the programs and activate the Conda environment that contains the program.
```bash
conda create -y -n seqQC -c conda-forge -c bioconda -c defaults cutadapt fastqc trimmomatic multiqc
```
Now activate your environment.
```bash
conda activate seqQC
```

Now to run the program we do the following, calling the sequence reads from the rawseq/16s directory like so but having the output placed in our current working directory:
```bash
fastqc ../../data_dir/*.fq -o .
```

The output is an .html file that you can view. You can look at them individually like so if mac user:
```bash
open B1_sub_R1_fastqc.html
```
Or like so if on mobaXterm
```bash
explorer.exe B1_sub_R1_fastqc.html
```


However, instead of checking each file individually, we can instead use the tool `multiqc` which will aggregate all of our results together.
```bash
multiqc .
open multiqc_report.html
```
Or like so if on mobaXterm
```bash
explorer.exe multiqc_report.html
```

Now we need to trim the primers off the reads and low-quality bases, both from read 1 and read 2, to improve their quality.


## 🧪 Exercise 2: Removing primers with Cutadapt 

Let's go back one directory, into the working_dir, and create a new directory called cutadapt.
```bash
cd ..
mkdir cutadapt
cd cutadapt/
```

Now copy over the data files we are working with to our current directory:
```bash
cp ../rawreads/*.fq .
```

Now we will run cutadapt on paired-end mode, because, remember from lecture, most cases reads are sequenced in the forward and reverse direction, meaning each forward read should have a paired read. We will try with one sample first before running as a loop to process all samples at once. 

```bash
cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
         -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
         -m 215 -M 285 --discard-untrimmed \
         -o B1_sub_R1_trimmed.fq -p B1_sub_R2_trimmed.fq \
         ../../data_dir/rawseqs/16s/B1_sub_R1.fq ../../data_dir/rawseqs/16s/B1_sub_R2.fq 
```

We are specifying the primers for the forward read with the -a flag, giving it the forward primer (in normal orientation), followed by three dots (required by cutadapt to know they are “linked”, with bases in between them, rather than right next to each other), then the reverse complement of the reverse primer. 

Then for the reverse reads, specified with the -A flag, we give it the reverse primer (in normal 5’-3’ orientation), three dots, and then the reverse complement of the forward primer. Both of those have a ^ symbol in front at the 5’ end indicating they should be found at the start of the reads (which is the case with this particular setup). The minimum read length (set with -m) and max (set with -M) were based roughly on 10% smaller and bigger than would be expected after trimming the primers. --discard-untrimmed states to throw away reads that don’t have these primers in them in the expected locations. Then -o specifies the output of the forwards reads, -p specifies the output of the reverse reads, and the input forward and reverse are provided as positional arguments in that order.

>[!NOTE]
> These types of settings will be different for data generated with different sequencing, i.e. not 2x300, and different primers sets. 

Ok, let's take a quick look to see that the primers were trimmed off. 
```bash
### R1 BEFORE TRIMMING PRIMERS
head -n 2 ../../data_dir/rawseqs/16s/B1_sub_R1.fq
# @M02542:42:000000000-ABVHU:1:1101:8823:2303 1:N:0:3
# GTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTCTTGT
# AAGACAGAGGTGAAATCCCTGGGCTCAACCTAGGAATGGCCTTTGTGACTGCAAGGCTGGAGTGCGGCAGAGGGGGATGG
# AATTCCGCGTGTAGCAGTGAAATGCGTAGATATGCGGAGGAACACCGATGGCGAAGGCAGTCCCCTGGGCCTGCACTGAC
# GCTCATGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCGGGTAGTCC

### R1 AFTER TRIMMING PRIMERS
head -n 2 B1_sub_R1_trimmed.fq
# @M02542:42:000000000-ABVHU:1:1101:8823:2303 1:N:0:3
# TACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTCTTGTAAGACAGAGGTGAAATCCC
# TGGGCTCAACCTAGGAATGGCCTTTGTGACTGCAAGGCTGGAGTGCGGCAGAGGGGGATGGAATTCCGCGTGTAGCAGTG
# AAATGCGTAGATATGCGGAGGAACACCGATGGCGAAGGCAGTCCCCTGGGCCTGCACTGACGCTCATGCACGAAAGCGTG
# GGGAGCAAACAGG


### R2 BEFORE TRIMMING PRIMERS
head -n 2 ../../data_dir/rawseqs/16s/B1_sub_R2.fq
# @M02542:42:000000000-ABVHU:1:1101:8823:2303 2:N:0:3
# GGACTACCCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCATGAGCGTCAGTGCAGGCCCAGGGGACTGCCTTC
# GCCATCGGTGTTCCTCCGCATATCTACGCATTTCACTGCTACACGCGGAATTCCATCCCCCTCTGCCGCACTCCAGCCTT
# GCAGTCACAAAGGCCATTCCTAGGTTGAGCCCAGGGATTTCACCTCTGTCTTACAAGACCGCCTGCGCACGCTTTACGCC
# CAGTAATTCCGATTAACGCTCGCACCCTACGTATTACCGCGGCTGCTGGCACTCACACTC

### R2 AFTER TRIMMING PRIMERS
head -n 2 B1_sub_R2_trimmed.fq
# @M02542:42:000000000-ABVHU:1:1101:8823:2303 2:N:0:3
# CCTGTTTGCTCCCCACGCTTTCGTGCATGAGCGTCAGTGCAGGCCCAGGGGACTGCCTTCGCCATCGGTGTTCCTCCGCA
# TATCTACGCATTTCACTGCTACACGCGGAATTCCATCCCCCTCTGCCGCACTCCAGCCTTGCAGTCACAAAGGCCATTCC
# TAGGTTGAGCCCAGGGATTTCACCTCTGTCTTACAAGACCGCCTGCGCACGCTTTACGCCCAGTAATTCCGATTAACGCT
# CGCACCCTACGTA
```

We are going to trim those again in the loop, so let's delete these trimmed files so as not to have duplicates.
```bash
rm *.fq
ls
```

Now, on to doing them all with a loop, here is how we can run it on all our samples at once. Since we have a lot of samples here, I’m redirecting the “stdout” (what’s printing the stats for each sample) to a file called *cutadapt_primer_trimming_stats.txt* so we can more easily view and keep track of if we’re losing a ton of sequences or not by having that information stored somewhere – instead of just plastered to the terminal window. We’re also going to take advantage of another convenience of cutadapt – by adding the extension .gz to the output file names, it will compress the files for us.

```bash
nano cutadapt.sh
```

Add this to the bash file you just created.
```bash
#!/usr/bin/env bash

DIR="../../data_dir/rawseqs/16s"

for R1 in ${DIR}/*_R1.fq; do
    # derive sample name by stripping path and suffix
    sample=$(basename "$R1" _R1.fq)
    R2="${DIR}/${sample}_R2.fq"

    echo "On sample: $sample"

    cutadapt \
        -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
        -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
        -m 215 -M 285 --discard-untrimmed \
        -o ${sample}_R1_trimmed.fq.gz \
        -p ${sample}_R2_trimmed.fq.gz \
        "$R1" "$R2" \
        >> cutadapt_primer_trimming_stats.txt 2>&1
done
```

- a = forward adapter
- A = reverse adapter
- m 215 = will discard all reads below 215 bp 
- M 285 = will discard all reads larger than 285 bp
- discard-untrimmed = discards reads in which no adapter match was found.
- o = output R1 
- p = output R2

Now run it
```bash
bash cutadapt.sh
ls
```

You should see all trimmed files here and the output stats file.
I typically like to have a file with all the sample names to use for various things throughout, so here’s making that file based on how these sample names are formatted. 
```bash
ls *_R1_trimmed.fq.gz | cut -f1 -d "_" > samples.txt
```

Now you can look through the output of the cutadapt stats file we made (“cutadapt_primer_trimming_stats.txt”) to get an idea of how things went. Here’s a little one-liner to look at what fraction of reads were retained in each sample (column 2) and what fraction of bps were retained in each sample (column 3):
```bash
paste samples.txt <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")
```

Great, so in all cases >90% of reads were kept. 

Let's also check again with FastQC/MultiQC to see how that improved the output
```bash
cd ../fastqc/
mkdir trimmed
cd trimmed/
fastqc ../../cutadapt/*.fq.gz -o .
multiqc .
open multiqc_report.html
```
With primers removed, we’re now ready to switch to R Studio and start using DADA2!

But, before we do that, let's prepare a new directory to work in with everything we need.
```bash
cd ../../
mkdir dada2
cd dada2
cp ../cutadapt/samples.txt .
ln -s ../cutadapt/*.fq.gz .
ls
```
I am introducing a new command here `ln`. This command creates a symbolic (hence -s) link, also known as a symlink or soft link. This is a special type of file that points to another file or directory. Symbolic links are commonly used to create shortcuts or aliases for files or directories located in the file system. This allows us to use those files in this directory, without having to make a duplicate, hard copy. 

Ok, now we can switch to R and process these reads (: 
