# Week 5: Gene-centric metagenomics

In this tutorial, we will annotate and analyze the genes from metagenome samples of the infant gut, published by Sharon et al 2023 (doi: 10.1101/gr.142315.112). 
The goal was to track species and strain level variations in microbial communities across 11 fecal samples collected from a premature infant during the first month of life. 
We are adapting a tutorial made for Anvio by Murat Eren (thank you 😊) 

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
- Understand the difference between a genome and a metagenome.
- Compare genes across metagenome samples. 

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

> If you type `ls`in the dataset directory you will see that the data-pack contains an anvi’o contigs database, an anvi’o merged profile database (that describes 11 metagenomes), and other additional data that are required by various sections in this tutorial. Here are some simple descriptions for some of these files, and how they were generated.
> An anvi’o contigs database was generate using the program `anvi-gen-contigs-database`. This special anvi’o database keeps all the information related to your contigs: positions of open reading frames, k-mer frequencies for each contig, functional and taxonomic annotation of genes, etc. The contigs database is an essential component of everything related to anvi’o metagenomic workflow.
> A merged anvi’o profile database was also generated using the program `anvi-profile`. In contrast to the contigs database, anvi’o profile databases store sample-specific information about contigs. Profiling a BAM file with anvi’o creates a single profile that reports properties for each contig in a single sample based on mapping results.
> Each profile database automatically links to a contigs database, and anvi’o can merge single profiles that link to the same contigs database into an anvi’o merged profile (which is what you will work with during this tutorial), using the program `anvi-merge`.
>
> Then the program `anvi-run-hmms` was used to identify single-copy core genes (SCGs) for Bacteria, Archaea, Eukarya, as well as sequences for ribosomal RNAs among the contigs. All of these results are also stored in the contigs database. This information allows us to learn the completion and redundancy estimates of newly identified bins (to be done next week) in the interactive interface, on the fly.
> To assign functions, the `anvi-run-ncbi-cogs` and `anvi-run-kegg-kofams` was run on the contigs database before, which stored functions for genes results in the contigs database. At the end of the binning process, functions occurring in each bin will be made available for downstream analyses.



Back to the tutorial, navigate to the `work_dir` and activate the anvio conda environment. 
```bash
cd ../../work_dir/
conda activate anvio-8
```

Make a shortcut (symlink) to the files here.
```bash
ln -s ../data_dir/INFANT-GUT-TUTORIAL/* .
ls
```

## 🧪 Step 2: Data analysis

Let’s take a first look at the merged profile database for the infant gut dataset metagenome. If you want to know more or have a reference to the Anvio interface, then follow this link: https://merenlab.org/tutorials/interactive-interface/
```bash
anvi-interactive -p PROFILE.db -c CONTIGS.db
```
<img width="2000" height="1730" alt="image" src="https://github.com/user-attachments/assets/81606879-9c2f-44d7-bdbc-61a6ce68c650" />



Important taxonomy 
```bash
anvi-import-taxonomy-for-genes -c CONTIGS.db \
                               -p centrifuge \
                               -i additional-files/centrifuge-files/centrifuge_report.tsv \
                               additional-files/centrifuge-files/centrifuge_hits.tsv
```
```
anvi-interactive -p PROFILE.db -c CONTIGS.db
```
Once gene-level taxonomy is added into the contigs database, anvi’o will determine the taxonomy of each contig based on the taxonomic affiliation of genes they describe, and display them in the interface whenever possible.

So at this point we don’t have any idea about what genomes do we have in this dataset, but anvi’o can make sense of the taxonomic make up of a given metagenome by characterizing taxonomic affiliations of single-copy core genes. 

```bash
anvi-estimate-scg-taxonomy -c CONTIGS.db --metagenome-mode
```

Nice but not as informative as it could be. Let's take use of our profile database the following way, which will give us a little more information about our dataset (remember its a time series:
```bash
anvi-estimate-scg-taxonomy -c CONTIGS.db \
                           -p PROFILE.db \
                           --metagenome-mode \
                           --compute-scg-coverages
```

```bash
$ anvi-display-contigs-stats contigs.db
```

```bash
anvi-interactive -c CONTIGS.db -p PROFILE.db
```

