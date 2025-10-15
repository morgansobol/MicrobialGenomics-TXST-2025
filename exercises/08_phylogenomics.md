# Week 7: Phylogenomics

In this tutorial, we will...

Putting genomes in a phylogenomic context is one of the common ways to compare them to each other. 
The common practice is to concatenate aligned sequences of single-copy core genes for each genome of interest, and generate a phylogenomic tree by analyzing the resulting alignment.

## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
* Explain


## 🧪 Step 1: Reading in the data and setting up the working environment

Like the last few weeks, you need to make a new `week8` directory inside your `microgenomics-2025` directory. 
Download the `work_dir` for this week from Canvas under Module Week 8 into 'week8` directory. 
To unpack it:
```
unzip work_dir.zip
```
Move into `work_dir`.

Last thing, don't forget to activate your anvio conda environment. 

## 🧪 Step 2: Ensure you have the correct bins
```
anvi-show-collections-and-bins -p PROFILE.db

```
You should see at least the "default" collection:
```
Collection: "default"
===============================================
Collection ID ................................: default
Number of bins ...............................: 13
Number of splits described ...................: 4,451
Bin names ....................................: Aneorococcus_sp, C_albicans, E_facealis, F_magna, L_citreum, P_acnes, P_avidum,
                                                P_rhinitidis, S_aureus, S_epidermidis, S_hominis, S_lugdunensis, Streptococcus
```

If not, run this to import it:
```
anvi-import-collection additional-files/collections/merens.txt \
                       --bins-info additional-files/collections/merens-info.txt \
                       -p PROFILE.db \
                       -c CONTIGS.db \
                       -C default
```

## 🧪 Step 3: Select genes from an HMM profile
In order to do the phylogenomic analysis, we will need a FASTA file of concatenated SCGs. And to get that FASTA file out of our anvi’o databases, we will primarily use the program `anvi-get-sequences-for-hmm-hits`.
We first need to identify an HMM profile, and then select some gene names from this profile to play with. Let’s start by looking at what HMM profiles are available to us: 
```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -p PROFILE.db \
                                -o seqs-for-phylogenomics.fa \
                                --list-hmm-sources
```
`Bacteria_71` is a collection that anvi’o developers curated by taking Mike Lee’s bacterial SCGs collection, first released in GToTree, which is another easy-to-use phylogenomics workflow.
Let’s see what genes do we have in `Bacteria_71`:
```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -p PROFILE.db \
                                -o seqs-for-phylogenomics.fa \
                                --hmm-source Bacteria_71 \
                                --list-available-gene-names
```

For the sake of this simple example, let’s assume we want to use a bunch of ribosomal genes for our phylogenomic analysis: Ribosomal_L1, Ribosomal_L2, Ribosomal_L3, Ribosomal_L4, Ribosomal_L5, Ribosomal_L6.
The following command will give us all these genes from all bins described in the collection `default`:
```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -p PROFILE.db \
                                -o seqs-for-phylogenomics.fa \
                                --hmm-source Bacteria_71 \
                                -C default \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6
```

Let's view the output:
```
less seqs-for-phylogenomics.fa
```
Press `q` to exit out of `less`. 


You can see this is not an alignment file but we need a single concatenated alignment of gene sequences per genome. 

Actually, the same function `anvi-get-sequences-for-hmm-hits` has a flag, `--concatenate-genes` to get your genes of interest to be concatenated. 
Let’s do that:
```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -p PROFILE.db \
                                -o seqs-for-phylogenomics.fa \
                                --hmm-source Bacteria_71 \
                                -C default \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
                                --concatenate-genes

```

Now, you are going to get an error that we need to use the flag `--return-best-hit`. The reason being that even in the most complete genomes, there may be multiple HMM hits for a given SCG.
As a solution to this problem, anvi’o asks you to use the --return-best-hit flag, which will return the most significant HMM hit if there are more than one gene that matches to the HMM of a given gene in a given genome. Fine. 
Let’s do that, then:
```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -p PROFILE.db \
                                -o seqs-for-phylogenomics.fa \
                                --hmm-source Bacteria_71 \
                                -C default \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
                                --concatenate-genes \
                                --return-best-hit
```

Let's look at the output again
```
less seqs-for-phylogenomics.fa
```

We are getting somewhere, but not quite there yet. See, the output is in DNA alphabet, which may not be the best option for phylogenomic analyses, especially if the genomes you have are coming from distant clades (which happens to be the case for IGD). 
Fortunately, you can easily switch to AA alphabet with an additional flag `--get-aa-sequences`. Ok, one last time:
```
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
                                -p PROFILE.db \
                                -o seqs-for-phylogenomics.fa \
                                --hmm-source Bacteria_71 \
                                -C default \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
                                --concatenate-genes \
                                --return-best-hit \
                                --get-aa-sequences
```
```
less seqs-for-phylogenomics.fa
```
Congrats, we did it. Moving on. 

## 🧪 Step 4: Compute a phylogenomic tree

Once you have your concatenated genes, which you now have them in `seqs-for-phylogenomics.fa`, it is time to perform the phylogenomic analysis.
Here we will use the program `anvi-gen-phylogenomic-tree`, which accepts a FASTA file and uses FastTree to infer phylogenetic trees.
```
anvi-gen-phylogenomic-tree -f seqs-for-phylogenomics.fa \
                           -o phylogenomic-tree.txt
```

Visualize the tree like so:
```
anvi-interactive --tree phylogenomic-tree.txt \
                 -p temp-profile.db \
                 --title "Pylogenomics of IGD Bins" \
                 --manual
```

Now, we can see how related our bins are. For instance, we can see that the S. hominis, S. epidermidis, and S. aureus, all Staphylococci, are correctly placed next to one another. 
The remaining bins are different organisms, placed appropriately from one another. 

## 🧪 Step 5: Phylogenomics for closely related organisms

What makes ribosomal SCGs powerful for some applications of phylogenomics, makes them weak for other applications. Comparing very similar genomes (pangenomics) may require the inclusion of a larger number of genes that occur in all genomes of interest to increase the signal of divergence. 

Let's test this with our_ E. faecalis_ and _E. faecium_ genomes from last week. We will first grab the same ribosomal SCGs like before:

```
anvi-get-sequences-for-hmm-hits --external-genomes additional-files/pangenomics/external-genomes.txt \
                                -o concatenated-proteins.fa \
                                --hmm-source Bacteria_71 \
                                --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
                                --concatenate-genes \
                                --return-best-hit \
                                --get-aa-sequences
```

Compute the tree and visualize it
```
anvi-gen-phylogenomic-tree -f concatenated-proteins.fa \
                           -o phylogenomic-tree.txt
```
```
anvi-interactive -p phylogenomic-profile.db \
                 -t phylogenomic-tree.txt \
                 --title "Enterococcus genomes" \
                 --manual
```

As you can see, many of the genomes have a flat line in the phylogenomic tree, which indicates that they have an identical set of the ribosomal proteins we used. This could often be the case for very closely related populations you want to study. 

Let's use the pangenomic analysis from last week. Select 100 genes randomly, and store them in a new collection called "Some_SCGs" to generate a phylogenomic tree using those:
```
anvi-display-pan -g Enterococcus-GENOMES.db \
                 -p PAN/Enterococcus-PAN.db \
                 --title "Enterococccus Pan"
```

The program anvi-get-sequences-for-gene-clusters is what we will use to export alignments genes in gene clusters. Here, we declare the collection name and the bin id in our anvi’o pan database, and export sequences:
```
anvi-get-sequences-for-gene-clusters -g Salmonella-GENOMES.db \
                                     -p Salmonella/Salmonella-PAN.db \
                                     --collection-name default \
                                     --bin-id Some_Core_PCs \
                                     --concatenate-gene-clusters \
                                     -o concatenated-proteins.fa
```




