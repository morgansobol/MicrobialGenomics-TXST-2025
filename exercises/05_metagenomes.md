# Week 5: Gene-centric metagenomics

In this tutorial, we will explore the presence of _Enterococcus faecalis_ across metagenome samples from the Human Microbiome Project (HMP), specifically infant guts (pub here: https://www.nature.com/articles/s41467-017-02018-w). 

_E. faecalis_ is a Gram-positive, commensal bacterium naturally inhabiting the gastrointestinal tracts of humans. Like other species in the genus Enterococcus, E. faecalis is typically found in healthy humans but it is also an opportunistic pathogen capable of causing severe infections, especially in a hospital setting.(https://en.wikipedia.org/wiki/Enterococcus_faecalis) 
We are adapting a tutorial made for Anvio by Murat Eren (thank you 😊) 

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
- Understand the difference between a genome and a metagenome.
- Assess differential coverage and detection of genes across metagenomes. 
- Compare genes across metagenome samples.
- Assesses differences in SNVs across genes.

## 🧪 Step 1: Reading in the data and setting up the working environment

1. Navigate to the `microgenomics-2025` directory you made last week.
2. Make a new directory called `week5` and move into that directory.
3. Now make two new directories, `data_dir` and `work_dir`.
4. Move into the `data_dir`.
5. Download and unpack the data like so:

```bash
curl -L -o INFANT-GUT-TUTORIAL.tar.gz \
     -H "User-Agent: Chrome/115.0.0.0" \
     https://figshare.com/ndownloader/files/45076909
```
```bash
tar -zxvf INFANT-GUT-TUTORIAL.tar.gz && cd INFANT-GUT-TUTORIAL
```

If you type `ls`in the dataset directory, you will see that the data-pack contains a lot of data. A lot of this we will use next week, but let's go ahead navigate to the `work_dir` to get the data we need today and activate the anvio conda environment. 
```bash
cd ../../work_dir/
conda activate anvio-8
```

Make a shortcut (symlink) to the files we will use today here.
```bash
ln -s ../data_dir/INFANT-GUT-TUTORIAL/additional-files/e_faeealis_across_hmp/* .
ls
```

> You should see an anvi’o contigs database of our reference genome, an anvi’o merged profile database that describes 20 gut metagenomes, and other additional data that are required by various sections in this tutorial. Here are some simple descriptions for some of these files and how they were generated.
> 
> An anvi’o contigs database was generated using the program `anvi-gen-contigs-database`. This special anvi’o database keeps all the information related to your contigs: positions of open reading frames, k-mer frequencies for each contig, functional and taxonomic annotation of genes, etc. The contigs database is an essential component of everything related to anvi’o metagenomic workflow.
> 
> A merged anvi’o profile database was also generated using the program `anvi-profile`. In contrast to the contigs database, anvi’o profile databases store sample-specific information about contigs. A mapping/alignment tool (like Bowtie) was used to map the metagenome short reads to our reference, creating a BAM file. Profiling a BAM file with anvi’o creates a single profile that reports properties for each contig in a single sample based on mapping results.
> Each profile database automatically links to a contigs database, and anvi’o can merge single profiles that link to the same contigs database into an anvi’o merged profile (which is what you will work with during this tutorial), using the program `anvi-merge`.
>
> When generating a profile with `anvi-profile`, we can ask it to profile the SNVs in our reads.
> 
> To assign functions, the `anvi-run-ncbi-cogs` was run on the contigs database before, which stored functions for genes in the contigs database. 


## 🧪 Step 2: Gene presence-absence

First, this data was prepared with an older version of Anvi'o, so we need to migrate the data to work with our version like so:
```bash
anvi-migrate --migrate-safely *.db
```

Let’s take a first look at the dataset and put our _E. faecalis_ contigs in the context of 20 infant gut metagenomes.
```bash
anvi-interactive -p PROFILE.db -c CONTIGS.db
```

We have some additional data to add that tells us more about the sample:
```bash
anvi-import-misc-data additional_layers_table.txt -p PROFILE.db -t layers
```
```bash
anvi-interactive -p PROFILE.db -c CONTIGS.db
```

Instead of doing this per-contig basis, we can look at the distribution patterns of genes. Don't worry, you will see warnings because it's our first time running this, and it needs to generate the gene database for us.
```bash
anvi-interactive -p PROFILE.db -c CONTIGS.db -C DEFAULT -b EVERYTHING --gene-mode
```
If you right-click any of the genes, you can see that this menu has more options, including viewing gene context. 


Let's output some data files on gene coverage and detection for you to use for the assignment:
```bash
anvi-export-gene-coverage-and-detection -c CONTIGS.db -p PROFILE.db -O efae
```

## 🧪 Step 3: SNVs

Now let's assess SNVs. As you can imagine, there are a lot of SNVs across all this data. We can get some context by looking at particular genes. 
```bash
anvi-interactive -p PROFILE.db -c CONTIGS.db -C DEFAULT -b EVERYTHING --gene-mode
```

We can export the data, but it would be just too much. So instead let's look at SAAVS (perhaps a bit more meaningful). 
```
anvi-gen-variability-profile -c CONTIGS.db -p PROFILE.db -C DEFAULT -b EVERYTHING --engine NT --min-coverage-in-each-sample 1 --quince-mode --compute-gene-coverage-stats -o snvs.txt 
```

Some info on settings:
* --engine = will define the output profile you will get from this program. The engine can focus on nucleotides (NT), codons (CDN), or
  an amino acids (AA).
* --min-coverage-in-each-sample = Minimum coverage of a given variable nucleotide position in all samples. If a nucleotide position is covered less than this value even in one
                        sample, it will be removed from the analysis. Default is 0.
* --quince-mode = The default behavior is to report allele frequencies only at positions where variation was reported during profiling (which by default uses
                        some heuristics to minimize the impact of error-driven variation). So, if there are 10 samples, and a given position has been reported as a
                        variable site during profiling in only one of those samples, there will be no information will be stored in the database for the remaining 9.
                        When this flag is used, we go back to each sample, and report allele frequencies for each sample at this position, even if they do not vary.
                        It will take considerably longer to report when this flag is on, and the use of it will increase the file size dramatically, however it is
                        inevitable for some statistical approaches and visualizations.
* --compute-gene-coverage-stats = If provided, gene coverage statistics will be appended for each entry in variability report. This is very useful information, but will not be
                        included by default because it is an expensive operation, and may take some additional time.


Finally, so you can link where SAAVs are to function, output the assigned gene functions to a file and the fastas for each gene
```bash
anvi-export-functions -c CONTIGS.db -o cog_functions.txt --annotation-sources COG14_CATEGORY,COG14_FUNCTION
anvi-get-sequences-for-gene-calls -c contigs.db -o genes_nt.fa
```


## 📝 Assignment due next class on Canvas
1. Do you find genes that are present in C-section, not present in vaginal birth samples?
2. Which genes have the most AA variation? Normalize by gene length because longer genes have more sites, so raw counts SAAVs will look bigger just due to length, not biology. Normalizing by length lets you compare apples to apples.
3. What are these genes and what functions are involved in? (Use fasta and cog_functions.txt)
