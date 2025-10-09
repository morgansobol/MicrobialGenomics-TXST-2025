# Week 7: Pangenomics

In this tutorial, we will be using a different part of the Infant Gut Dataset containing 11 feacal samples from a premature infant collected between Day 15 and Day 24 of her life in a time-series manner. 

## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
* Explain the concept of a pangenome and distinguish between core, accessory, and unique genes across a set of genomes.
* Use Anvi’o to generate and visualize a pangenome.
* Interpret patterns in sequence similairty and gene presence/absence matrices.


## 🧪 Step 1: Reading in the data and setting up the working environment

Like the last few weeks, you need to make a new `week7` directory, with a `work_dir` inside your `microgenomics-2025` directory. 

## 🧪 Step 2: Initial visualization of pangenome

For this example, we will work with 6 external _E. faecalis_, and 5 _E. faecium_ genomes to analyze them together with our _E. faecalis_ bin from last week. For each of these 11 external genomes, anvi’o contigs databases were already created by Meren. 

We can see them here:
```
ls additional-files/pangenomics/external-genomes/*db
```

The first step in the pangenomic workflow is to generate an anvi’o genomes storage, which is a special anvi’o database that stores information about genomes. Easy enough!
```
anvi-gen-genomes-storage -i additional-files/pangenomics/internal-genomes.txt \
                         -e additional-files/pangenomics/external-genomes.txt \
                         -o Enterococcus-GENOMES.db
```

Now, we can characterize the pangenome. Conveniently, Anvi'o has a `pan-genome` command for us. 
```
anvi-pan-genome -g Enterococcus-GENOMES.db \
                -n Enterococcus \
                -o PAN \
                --num-threads 4
```

Now the pangenome is ready to display. This is how you can run it:
```
anvi-display-pan -g Enterococcus-GENOMES.db \
                 -p PAN/Enterococcus-PAN.db \
                 --title "Enterococccus Pan"
```

Meren made the visualization better for us; let's import his "state" file. Kill the server, and import the state file the following way, and re-run the server,
```
anvi-import-state -p PAN/Enterococcus-PAN.db \
                  --state additional-files/state-files/state-pan.json \
                  --name default
```
```
anvi-display-pan -g Enterococcus-GENOMES.db \
                 -p PAN/Enterococcus-PAN.db \
                 --title "Enterococccus Pan"
```

Now not only can we see how our _E. faecalis_SHARON MAG looks like compared to available genomes, but we can also see that it is not missing or carrying a great number of proteins compared to other genomes. The clustering of genomes based on gene clusters indicates that it is most similar to the genome of _Enterococcus faecalis_ 6512. 

## 🧪 Step 3: Average nucleotide identity 

The pangenome tells us about the similarities and dissimilarities between those genomes given the amino acid sequences of open reading frames we identified within each one of them. But we could also compare genomes to each other by computing the average nucleotide identity between them. Anvi’o comes with a program to do that, which uses PyANI by Pritchard et al. to pairwise align DNA sequences between genomes to estimate similarities between them. 

We can compute average nucleotide identity among our genomes, and add them to the pan database in the following way:
```
anvi-compute-genome-similarity -e additional-files/pangenomics/external-genomes.txt \
                               -i additional-files/pangenomics/internal-genomes.txt \
                               --program pyANI \
                               -o ANI \
                               -T 4 \
                               --pan-db PAN/Enterococcus-PAN.db
```
This can take a bit of time, especially as you add more genomes. 

Once it's done, visualize it again
```
anvi-display-pan -g Enterococcus-GENOMES.db \
                 -p PAN/Enterococcus-PAN.db \
                 --title "Enterococccus Pan"
```

And take a careful look at your ‘Layers’ tab, you can visualize your pangenome with the new ANI results.

## 🧪 Step 4: Organizing gene clusters based on a forced synteny
By default, anvi’o will organize gene clusters in a given pangenome based on their distribution patterns across genomes. Which, by definition, will break gene synteny, i.e., two gene clusters that are next to each other in the display will not necessarily carry genes that are next to each other in respective genomes. If you have a relatively complete genome, you can ask anvi’o to organize gene clusters based on the organization of genes in that genome.

Here, for instance, we can force the organization of gene clusters based on the synteny of genes in our _E. faecalis_ genome:


## 🧪 Step 5: Functional enrichment analyses
One of the things we often are interested in is this question: which functions are associated with a particular organization of genomes due to their phylogenomic or pangenomic characteristics.

For instance, genomes in this pangenome are organized into two distinct groups based on differential distribution of gene clusters. That is of course not surprising, since these genomes are classified into two distinct ‘species’ within the genus Enterococcus. So you can certainly imagine more appropriate or interesting examples where you may be wondering about functional enrichment across groups of genomes that do not have such clear distinctions at the level of taxonomy.

This is done with `anvi-compute-functional-enrichment-in-pan`.

This program utilizes information you provide by splitting your genomes into multiple grups, and finds functions that are characteristic to any of those groups, i.e., functions that are enriched in genomes that belong to one group, but predominantly absent in genomes outside of it.

uses a Generalized Linear Model with the logit linkage function to compute an enrichment score and p-value for each function. False Detection Rate correction to p-values to account for multiple tests is done using the package qvalue.

In addition to the enrichment test, we use a simple heuristic to find the groups that associate with each function. This association is only meaningful for functions that are truly enriched, and should otherwise be ignored. We simply determine that for every function, the associated groups are the ones in which the occurrence of the function of genomes is greater than the expected occurrence under a uniformal distribution (i.e. if the function was equally probable to occur in genomes from all groups). 

First, you need to define how you would like to group your genomes in the pangenome so anvi’o can find out which functions are characteristic of each group. We'll group by _faecalis_ and _faecium_.

Once it is imported, we can re-run the pangenome,

anvi-display-pan -g Enterococcus-GENOMES.db \
                 -p PAN/Enterococcus-PAN.db \
                 --title "Enterococccus Pan"

                 Now we can use the program anvi-compute-functional-enrichment-in-pan to ask anvi’o to identify and report functions that are enriched in either of these clades along with the gene clusters they are associated with:

anvi-compute-functional-enrichment-in-pan -p PAN/Enterococcus-PAN.db \
                                          -g Enterococcus-GENOMES.db \
                                          --category-variable clade \
                                          --annotation-source COG20_FUNCTION \
                                          -o functional-enrichment.txt


