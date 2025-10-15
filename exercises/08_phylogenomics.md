# Week 7: Phylogenomics

In this tutorial, we will...

Putting genomes in a phylogenomic context is one of the common ways to compare them to each other. 
The common practice is to concatenate aligned sequences of single-copy core genes for each genome of interest, and generate a phylogenomic tree by analyzing the resulting alignment.

## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
* Explain


## 🧪 Step 1: Reading in the data and setting up the working environment

Like the last few weeks, you need to make a new `week8` directory inside your `microgenomics-2025` directory. 
Move into your `week8` directory and make a `work_dir`. Copy the contents of `week7/work_dir/` here.
```
ln -s ln -s ../../week7/work_dir/* .
```

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

## 🧪 Step 4: Compute a phylogenomic tree.

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

