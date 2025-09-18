# Week 4: Isolate/Whole genome assembly & annotation

In this tutorial, we are going to assemble an unknown genome and analyze it to find out what it is and what it does.
We will be back working in the command line, not in R.


---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
- 

## 🧪 Step 1: Setting up the working environment and reading in the data

Let's make a new `conda` environment for assembly. Open `Terminal` and run. 

```bash
conda create -n assembly -c bioconda -c conda-forge spades seqkit
```
If `conda` is not installed: you need to install following these instructions: https://github.com/morgansobol/MicrobialGenomics-TXST-2025/blob/main/exercises/install_conda.md 


Now let's make an environment for Anvio. 

1. Open this link: https://anvio.org/install/macos/stable/
2. Start at Step 3. Yes, your Mac has an Apple Silicon processor, so you need to run `conda config --env --set subdir osx-64` before creating the conda environment.
3. Stop at Step 4. If you receive no errors, proceed to Step 6 and run `anvi-self-test --suite mini` to check your installation.

Ok, hopefully that all worked. We will now create our directories and download the data. 

Let's start with a fresh directory. Make sure you are in your ~home directory, run `cd` if not.

Run these lines one at a time. 
```bash
cd Desktop
mkdir microgenomics-2025
cd microgenomics-2025
mkdir week4
cd week4
mkdir data_dir working_dir
cd data_dir
```

Now download the two fastq files, `unknown_R1_paired.fastq.gz` and `unknown_R2_paired.fastq.gz` from Canvas to Downloads and move them over to `data_dir`. 
```bash
mv ../../../../Downloads/*.fastq.gz .
ls
```

Now let's move ourselves into `work_dir` and start processing our reads there. Don't forget to activate the conda environment.
```bash
cd ../work_dir/
```

>[!NOTE]
> We are going to skip read QC and trimming for the sake of time. So here is what the authors did.

Sequencing: Illumina Nextseq, 2x150bp, with a 350bp insert. 

<img width="800" height="600" alt="image" src="https://github.com/user-attachments/assets/ac70010c-774b-41c6-a1d4-8662f8ab2357" />

The authors of this tutorial used a different program called `Trimmomatic`. It's like cutadapt, but a bit more automatic. It can automatically detect Illumina adapters.
This was the code they ran:
```bash
trimmomatic PE unknown_raw_R1.fastq.gz unknown_raw_R2.fastq.gz \
            unknown_R1_paired.fastq.gz unknown_R1_unpaired.fastq.gz \
            unknown_R2_paired.fastq.gz unknown_R2_unpaired.fastq.gz \
            CROP:140 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 \
            MINLEN:140 -threads 4
```
* PE = run in paired-end mode
* LEADING:10 = cut the bases off the start of the read if their phred quality score is below 10
* TRAILING:10 = does the same, but on the ends of the reads
* SLIDINGWINDOW:5:20 = starting at base 1, look at a window of 5 bps and if the average quality score drops before 20, truncate the read at that position and only keep up to that point
* MINLEN:151 = pretty stringent, any read that gets trimmed will be thrown out

<img width="800" height="600" alt="image" src="https://github.com/user-attachments/assets/95df4988-d3c8-4e98-aa31-8f61a9bcb354" />


## 🧪 Exercise 2: Assembly

Make a directory for SPAdes to work in. 

```bash
mkdir spades
cd spades
```
Let's create a symlink to our data from data_dir.

```bash
ln -s ../../data_dir/*.fastq.gz .
```
Now we are ready to assemble!
The manual for SPAdes can be found here: https://ablab.github.io/spades/
We can also view it on our command line, run `spades.py -h`. 

The steps are: 

<img width="265" height="355" alt="Screenshot 2025-09-18 at 10 42 34 AM" src="https://github.com/user-attachments/assets/9ad8209b-6707-499e-b627-594461cb6f39" />
<img width="479" height="267" alt="Screenshot 2025-09-18 at 10 43 29 AM" src="https://github.com/user-attachments/assets/f65ab785-ac47-4846-88a5-3ab6be2075f4" />



On their manual, SPAdes recommends a few things for us:
1. We are working with an isolate genome (not metagenome), so the read coverage across the genome is probably pretty high (>50x). We need to run our assembly in `--isolate` mode.
2. For 150bp reads, they recommend kmers of 21,33,55,77.

> Small k → better sensitivity, connects through low-coverage regions but introduces tangles/repeats.
> Large k → more specificity, helps resolve repeats, but risks breaking contigs in low-coverage areas.

~Rule of thumb: the largest k should be about ½ the read length or a bit more.~

You may ask, why the odd numbers? 
--> With an even k, you can end up with k-mers that are perfect palindromes, i.e. a sequence that reads the same in forward and reverse in the reverse complement of the sequence. 
```
                5' - TCGCGA - 3'
                3' - AGCGCT - 5'
                5' - TCGCGA - 3'

                5' - TCGCG - 3'
                3' - AGCGC - 5'
                5' - CGCGA - 3'
```
        
Ok, I think that covers SPAdes. We can finally assemble. It only requires one command from us:            
>[!WARNING] this step consumes at least 16 GB of RAM/memory. If you have an 8 GB computer, SPAdes will likely fail since you have other things running. 

```bash
spades.py --isolate -1 unknown_R1_paired.fastq.gz -2 unknown_R2_paired.fastq.gz -o output -t 4 -k 21,33,55,77
```
This will take <5 min. 

My favorite assembly stats program, Quast, was not working 😔
So we use seqkit's `stat` function instead. 

```bash
seqkit stats -N 50 output/contigs.fasta
```

We can also use a bash command `grep` to grab and then count (-c) how many contigs we have by using the `>`. 
```bash
grep -c ">" output/contigs.fasta
```

Now we’re going to put our genome assembly into the anvi’o framework and begin to look at our assembled genome.

Deactivate your conda environment.
```
conda deactivate
```

## 🧪 Exercise 3: Anvi'o

First, let's move back to our working_dir and make a new directory for anvi'o.
```bash
cd ../
pwd
mkdir anvio
cd anvio
```

Symlink the contigs.fasta and fastq files here:
```
ln -s ../spades/output/contigs.fasta .
ln -s ../../data_dir/*.fastq.gz .
```

Activate the conda environment.
```
conda activate anvio-8
```

For us to get our assembly into anvi’o, first we need to generate what it calls a contigs database using the `anvi-gen-contigs-database ` command. This will organize our contigs in an anvi’o-friendly way, and provide information about them. 

When run on `contigs.fasta` this program will:

1. Compute k-mer frequencies for each contig (the default is 4, but you can change it using --kmer-size parameter if you feel adventurous).
2. Soft-split contigs longer than 20,000 bp into smaller ones (you can change the split size using the --split-length flag). When the gene calling step is not skipped, the process of splitting contigs will consider where genes are and avoid cutting genes in the middle. For very, very large assemblies this process can take a while, and you can skip it with --skip-mindful-splitting flag.
3. Identify open reading frames using Prodigal, UNLESS, (1) you have used the flag --skip-gene-calling (no gene calls will be made) or (2) you have provided external-gene-calls.

You can see what the command needs and options you want to set by `anvi-gen-contigs-database -h `.

```bash
anvi-gen-contigs-database -f contigs.fasta -o contigs.db -n unknown_genome
```

```bash
anvi-run-hmms -c contigs.db -I Bacteria_71 -T 4
```
```bash
anvi-run-hmms -c contigs.db -I Ribosomal_RNA_16S -T 4
```

Anvio will determine our genome completeness and contamination with single-copy genes (SCGs) by running this:
```
anvi-estimate-genome-completeness -c contigs.db
```

So you can see, our genome is 100% complete based on the presence of expected SCGs. The redundancy is ~4%, so just below the 5% threshold for high-quality genomes. 

Now, let's assign functions to our ORFS!
First, we will do so with the NCBI COG database. We have to run `anvi-setup-ncbi-cogs` to download and set up the database. You only have to do this once. 

```bash
anvi-setup-ncbi-cogs -T 4 
anvi-run-ncbi-cogs -c contigs.db -T 4
```

Next, we will also assign function using the KEGG database. 
```bash
anvi-run-kegg-kofams -c contigs.db -T 4
```
This will take ~8-10 min.

When this next program runs, it will look at the KOfam annotations (for KEGG Orthologs, or KOs) within each genome, match them up to the KEGG module definitions to estimate the completeness of each module/pathway. 
A module is considered ‘complete’ or ‘present’ in a genome if its completeness score is above a certain threshold, which can be set with the --module-completion-threshold parameter. A static threshold such as this is not the most ideal metric, especially since metabolic modules have variable numbers of genes - for example, with the default threshold of 0.75 (75%), a module with 3 KOs in it would only be considered complete if all 3 of those KOs were found in a genome, while a module with 5 KOs could be considered complete if only 4 of its KOs were found. But it is what it is without diving deeper and doing additional checks (: 
```
anvi-estimate-metabolism -c contigs.db
```

We can then integrate mapping information from recruiting our reads to the assembly. This mapping information is placed into anvi’o with the anvi-profile program, which generates another type of database anvi’o calls a “profile database”. In contrast to the contigs-db, an anvi’o single-profile-db stores sample-specific information about contigs. Profiling a BAM file with anvi’o using anvi-profile creates a single profile that reports properties for each contig in a single sample based on mapping results. 

SAM files are a type of text file format that contains the alignment information of various sequences that are mapped against reference sequences. BAM files contain the same information as SAM files, except they are in binary file format which is not readable by humans. On the other hand, BAM files are smaller and more efficient for software to work with than SAM files, saving time and reducing costs of computation and storage. 

Run these line by line. 
```bash
bowtie2-build contigs.fasta contigs.btindex

bowtie2 -q -x contigs.btindex \
        -1 unknown_R1_paired.fastq.gz \
        -2 unknown_R2_paired.fastq.gz \
        -p 4 -S unknown_assembly.sam

samtools view -bS unknown_assembly.sam > unknown_assembly.bam

anvi-init-bam unknown_assembly.bam -o unknown_anvio.bam

anvi-profile -i unknown_anvio.bam -c contigs.db -T 4 --cluster-contigs -o unknown_profiled/

```

Let's export a fasta file of Prodigal-identified open-reading frames. 
```bash
anvi-get-sequences-for-gene-calls -c contigs.db -o gene_calls.fa
```

And do the same for our 16S rRNAs, SCGs, COGs, and KEGGs. Run each one line by line. 
```bash
anvi-get-sequences-for-hmm-hits -c contigs.db --hmm-source Ribosomal_RNA_16S -o rRNAs.fa
anvi-get-sequences-for-hmm-hits -c contigs.db --hmm-sources Bacteria_71 --get-aa-sequences -o bacterial_SCGs.faa --no-wrap
anvi-export-functions -c contigs.db -o cog_functions.txt --annotation-sources COG20_CATEGORY
anvi-export-functions -c contigs.db -o kegg_functions.txt --annotation-sources KEGG_Class,KOfam
```

Now we are going
```bash
anvi-script-add-default-collection -p unknown_profiled/PROFILE.db
anvi-summarize -c contigs.db -p unknown_profiled/PROFILE.db -C DEFAULT -o unknown_assembly_summary/
anvi-interactive -c contigs.db -p unknown_profiled/PROFILE.db --title "Unknown assembly"
```
Anvi'o visualization is really geared towards metagenomics/comparative genomics like so:
<img width="246" height="200" alt="image" src="https://github.com/user-attachments/assets/db450eb1-8b50-4d0d-b68e-8b791785cab5" />


## 📝 Assignment due next class on Canvas
1. Who does your genome belong to? (Hint: what can you do with the information in rRNAs.fa file)
2. Tell me about the organism (environments it's found in, metabolisms, etc). Provide literature references. 
3. Using the provided cog_functions.txt, make a bar chart of COG functional categories. Count the number of genes per COG category letter (A, C, E, …). If a gene is annotated to multiple categories (entries separated by !!!), count it once for each category it belongs to. Label axes and add category names. Feel free to use Excel (easy) or R (advanced). If you use ChatGPT, provide a screenshot of the solution you used. 
Example:
<img width="1000" height="500" alt="image" src="https://github.com/user-attachments/assets/8f0804d9-c542-45f6-9c4b-780c4f335de7" />









