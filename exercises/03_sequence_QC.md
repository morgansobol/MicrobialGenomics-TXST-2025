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

## 🧪 Exercise 1: FastQC, checking your read quality
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
fastqc ../../rawseqs/B1_sub_R1.fq ../../rawseqs/B1_sub_R2.fq
```
We should get an .html file to view the outputs.

[insert pics of output]

You can see the errors. Let's still trim both reads and see how they further improve. 


## 🧪 Exercise 2: Removing primers with Cutadapt 
Don't forget to go back one directory, into the working_dir, and create a new directory called cutadapt.
```bash
cd ..
mkdir cutadapt
cd cutadapt/
```
Now copy over the data files we are working with to our current directory:
```bash
cp ../rawreads/*.fq .
```

In our working directory there are 20 samples with forward (R1) and reverse (R2) reads with per-base-call quality information, so 40 fastq files (.fq). I typically like to have a file with all the sample names to use for various things throughout, so here’s making that file based on how these sample names are formatted. 
```bash
ls *_R1.fq | cut -f1 -d "_" > samples
```

Now we will run cutadapt on paired-end mode, because remember from lecture, most cases reads are sequenced in the forward and reverse direction, meaning each forward read should have a paired read. 
```bash
cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
         -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
         -m 215 -M 285 --discard-untrimmed \
         -o B1_sub_R1_trimmed.fq -p B1_sub_R2_trimmed.fq \
         B1_sub_R1.fq B1_sub_R2.fq 
```
We are specifying the primers for the forward read with the -a flag, giving it the forward primer (in normal orientation), followed by three dots (required by cutadapt to know they are “linked”, with bases in between them, rather than right next to each other), then the reverse complement of the reverse primer. 

Then for the reverse reads, specified with the -A flag, we give it the reverse primer (in normal 5’-3’ orientation), three dots, and then the reverse complement of the forward primer. Both of those have a ^ symbol in front at the 5’ end indicating they should be found at the start of the reads (which is the case with this particular setup). The minimum read length (set with -m) and max (set with -M) were based roughly on 10% smaller and bigger than would be expected after trimming the primers. --discard-untrimmed states to throw away reads that don’t have these primers in them in the expected locations. Then -o specifies the output of the forwards reads, -p specifies the output of the reverse reads, and the input forward and reverse are provided as positional arguments in that order.

>[!NOTE]
> These types of settings will be different for data generated with different sequencing, i.e. not 2x300, and different primers sets. 

Ok, let's take a quick look to see that the primers were trimmed off. 
```bash
### R1 BEFORE TRIMMING PRIMERS
head -n 2 B1_sub_R1.fq
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
head -n 2 B1_sub_R2.fq
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

Now, on to doing them all with a loop, here is how we can run it on all our samples at once. Since we have a lot of samples here, I’m redirecting the “stdout” (what’s printing the stats for each sample) to a file called *utadapt_primer_trimming_stats.txt* so we can more easily view and keep track of if we’re losing a ton of sequences or not by having that information stored somewhere – instead of just plastered to the terminal window. We’re also going to take advantage of another convenience of cutadapt – by adding the extension .gz to the output file names, it will compress them for us.

```bash
for sample in $(cat samples)
do

    echo "On sample: $sample"
    
    cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
             -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
             -m 215 -M 285 --discard-untrimmed \
             -o ${sample}_sub_R1_trimmed.fq.gz -p ${sample}_sub_R2_trimmed.fq.gz \
             ${sample}_sub_R1.fq ${sample}_sub_R2.fq \
             >> cutadapt_primer_trimming_stats.txt 2>&1

done
```

You can look through the output of the cutadapt stats file we made (“cutadapt_primer_trimming_stats.txt”) to get an idea of how things went. Here’s a little one-liner to look at what fraction of reads were retained in each sample (column 2) and what fraction of bps were retained in each sample (column 3):
```bash
paste samples <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")
```

Let's also check again with FastQC to see how that improved the output
```bash
cd ../fastqc/
fastqc ../cutadapt/R1.fq R2.fq
```
With primers removed, we’re now ready to switch to R and start using DADA2!
