# Week 6: Genome-centric metagenomics

In this tutorial, we will be using a different part of the Infant Gut Dataset containing 11 feacal samples from a premature infant collected between Day 15 and Day 24 of her life in a time-series manner. 

## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
- Familiarize yourself with the Anvi'o interface for binning. 
- Inspect contigs in the context of their metagenomic signal.
- Characterize bins by manual binning.
- Summarize manual binning results for downstream analyses.
- Manually curate individual bins for quality control. 
- Import and visualize external binning results.
- Combine manual and automatic binning.

## 🧪 Step 1: Reading in the data and setting up the working environment
1. Navigate to the `microgenomics-2025` directory.
2. Make a new directory called `week6` and move into that directory.
3. Let's just make a `work_dir` since we will use last week's `data_dir`. 
4. Move into the `work_dir` you just made.
5. Activate the `anvio-8` conda environment.
6. Make a shortcut (symlink) to the files we will use today here.
```bash
ln -s ../../data_dir/INFANT-GUT-TUTORIAL/* .
ls
```

## 🧪 Step 2: Inferring taxonomy
Let's take a first look at the merged profile database for the infant gut dataset metagenome. 
```bash
anvi-interactive -p PROFILE.db -c CONTIGS.db
```

Anvi’o can work with gene-level taxonomic annotations, but gene-level taxonomy is not useful for anything beyond occasional help with manual binning. Once gene-level taxonomy is added into the contigs database, anvi’o will determine the taxonomy of each contig based on the taxonomic affiliation of genes they describe, and display them in the interface whenever possible.

Import that data:
```bash
anvi-import-taxonomy-for-genes -c CONTIGS.db \
                               -p centrifuge \
                               -i additional-files/centrifuge-files/centrifuge_report.tsv \
                               additional-files/centrifuge-files/centrifuge_hits.tsv
```

```
anvi-interactive -p PROFILE.db -c CONTIGS.db
```
You will see an additional layer with taxonomy. In the Layers tab, find the Taxonomy layer, set its height to 200, then drag the layer in between DAY24 and hmms_Ribosomal_RNAs, and click Draw again. Then click the Save State button, and overwrite the default state. This will make sure anvi’o remembers to make the height of that layer 200px the next time you run the interactive interface!


At this point, we don’t have any idea about what genomes we have in this dataset, but anvi’o can make sense of the taxonomic make up of a given metagenome by characterizing taxonomic affiliations of single-copy core genes. 

Take a quick look at taxonomy:
```
anvi-estimate-scg-taxonomy -c CONTIGS.db --metagenome-mode
```

Good, but could be more informative. We can make use of our profile database in the following way, which will give us a little more information about our dataset:
```
anvi-estimate-scg-taxonomy -c CONTIGS.db \
                           -p PROFILE.db \
                           --metagenome-mode \
                           --compute-scg-coverages
```
Looks like information that would have been useful to have in front of us in our interactive interface. Luckily, anvi’o can add these taxonomic insights into a given profile database, if you change the previous command just a bit:
```
anvi-estimate-scg-taxonomy -c CONTIGS.db \
                           -p PROFILE.db \
                           --metagenome-mode \
                           --compute-scg-coverages \
                           --update-profile-db-with-taxonomy
```
```
anvi-interactive -c CONTIGS.db -p PROFILE.db
````
At this point, we have an overall idea about the makeup of this metagenome, but we don’t have any genomes from it. The following sections will cover some of the multiple ways to do this.

## 🧪 Step 3: Manual binning
I'm going to let you try to identify bins on your own first. A few tips:
* Contigs are already cluctered together based on tetranucleotide (k=4) frequency.
* Shared/similar coverage of contigs usually indicates contigs that belong together.
* If you click on the “Bins” tab at the top left, and then select the branch on the tree at the center that holds all the contigs, you will see a real-time estimate of % completion and redundancy.
* You do not have to bin all contigs. Instead, try to identify bins corresponding to an actual genome. 
* Please try to avoid bins with redundancy >10% to ensure they are more accurate and not contaminated.
* You can increase the inner tree radius (e.g., 5,000) for a better binning experience in the `Main` tab.
* You can select the option show grid in the Main tab (additional settings) for a better demarcation of identified bins.

Example:

<img width="2634" height="1184" alt="image" src="https://github.com/user-attachments/assets/391d9576-5d96-4053-9862-7d57f643f97f" />

Please save your bins as a `collection`. You can give your collection any name, but if you call it default, anvi’o will treat it differently. In the anvi’o lingo, a collection is something that describes one or more bins, each of which describes one or more contigs.

Let’s summarize the collection you have just created:
```
anvi-summarize -p PROFILE.db \
               -c CONTIGS.db \
               -C default \
               -o SUMMARY
```

As you can see from the summary file, at this point bin names are random. It would be more useful to put some order on this front. This becomes an extremely useful strategy, especially when the intention is to merge multiple binning efforts later. For this task we use the program anvi-rename-bins:
```
anvi-rename-bins -p PROFILE.db \
                 -c CONTIGS.db  \
                 --collection-to-read default  \
                 --collection-to-write MAGs  \
                 --call-MAGs \
                 --prefix IGD \
                 --report-file rename-bins-report.txt
```

With those settings, a new collection MAG will be created in which (1) bins with a completion >70% are identified as MAGs (stands for Metagenome-Assembled Genome), and (2) bins and MAGs are attached the prefix IGD and renamed based on the difference between completion and redundancy.

Now we can summarize the new collection:
```
anvi-summarize -p PROFILE.db \
               -c CONTIGS.db \
               -C MAGs \
               -o SUMMARY_AFTER_RENAME

```

You can now visualize the results by double-clicking on the index.html file present in the newly created folder SUMMARY_MAGs.

Refining individual MAGs




## 📝 Assignment due next class on Canvas
1. Calculate the relative abundance of bins across samples and generate a bar plot to show changes of bins across the days.
2. 
