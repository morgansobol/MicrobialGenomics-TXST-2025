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
4. `cd` into the `work_dir` you just made.
5. Activate the `anvio-8` conda environment.
6. Make a shortcut (symlink) to the files we will use today here.
```bash
ln -s ../../week5/data_dir/INFANT-GUT-TUTORIAL/* .
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
* Contigs are already clustered together based on tetranucleotide (k=4) frequency.
* Shared/similar coverage of contigs usually indicates contigs that belong together.
* You can increase the inner tree radius (e.g., 5,000) for a better binning experience in the `Main` tab, then save State and Draw. 
* You can select the option show grid in the Main tab (additional settings) for a better demarcation of identified bins.
* If you click on the “Bins” tab at the top left and then select the branch on the tree at the center that holds all the contigs, you will see a real-time estimate of % completion and redundancy.
* You do not have to bin all contigs. Instead, try to identify bins corresponding to an actual genome. 
* Please try to avoid bins with redundancy >10% to ensure they are more accurate and not contaminated.

Example:

<img width="2634" height="1184" alt="image" src="https://github.com/user-attachments/assets/391d9576-5d96-4053-9862-7d57f643f97f" />

Please save your bins as a `collection`. Let's give you collection a name "my_bins". In the anvi’o lingo, a collection is something that describes one or more bins, each of which describes one or more contigs.

Let’s summarize the collection you have just created:
```
anvi-summarize -p PROFILE.db \
               -c CONTIGS.db \
               -C my_bins \
               -o SUMMARY
```
```
open SUMMARY/index.html
```

As you can see from the summary file, at this point bin names are random. It would be more useful to put some order on this front. This becomes an extremely useful strategy, especially when the intention is to merge multiple binning efforts later. For this task we use the program anvi-rename-bins:
```
anvi-rename-bins -p PROFILE.db \
                 -c CONTIGS.db  \
                 --collection-to-read my_bins  \
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
```
open SUMMARY_AFTER_RENAME/index.html
```

## 🧪 Step 3: Refining our individual MAGs
To straighten the quality of the MAGs collection, it is possible to visualize individual bins and if needed, refine them. For this we use the program anvi-refine. For instance, if you were to be interested in refining one of the bins in our current collection, you could run this command:
```
anvi-refine -p PROFILE.db \
            -c CONTIGS.db \
            -C MAGs \
            -b IGD_MAG_00001
```
Now the interactive interface only displays contigs from a single bin. During this curation step, one can try different clustering strategies (i.e. by only relying on coverage, or only relying on sequence composition) to identify outliers and investigate carefully whether they may be contaminants. You can select everything and remove those contigs you don’t want to keep in the bin before using the Bins panel to store your updated set of contigs in the database.

Let's refine MAG_00001 together. I want you to do your other bins on your own. 
Storing the new refined bins in the database and it will modify the collection for you, but you will need to re-run `anvi-summarize` if you want the summary output to also be updated.

```
anvi-summarize -p PROFILE.db \
               -c CONTIGS.db \
               -C MAGs \
               -o SUMMARY_AFTER_REFINE
```

## 🧪 Step 3: Manual vs Automatic Binning
Even if you prefer manual binning over automatic binning for the sake of accuracy and control over your data, automatic binning is an unavoidable need due to performance limitations associated with manual binning

The directory additional-files/external-binning-results contains a number of files that describe the binning of contigs in the IGD based on various automatic and manual approaches. These files include (1) outputs from some of the well-known binning algorithms (i.e., GROOPM.txt, MAXBIN.txt, METABAT.txt, BINSANITY_REFINE.txt, MYCC.txt, and CONCOCT.txt), (2) the original binning of this dataset (SHARON_et_al.txt), and the manual binning performed in the anvi’o paper (MEREN_et_al.txt).

You can use the program anvi-show-collections-and-bins to see all collections in your anvi’o profile database:
```
anvi-show-collections-and-bins -p PROFILE.db
```
And use this one to remove one that’s called default, if you have it:
```
anvi-delete-collection -p PROFILE.db \
                       -C default
```
You can create a collection by using the interactive interface (e.g., the default and MAGs collections you just created), or you can import external binning results into your profile database as a collection and see how that collection groups contigs. For instance, let’s import the CONCOCT collection:
```
anvi-import-collection additional-files/external-binning-results/CONCOCT.txt \
                       -c CONTIGS.db \
                       -p PROFILE.db \
                       -C CONCOCT \
                       --contigs-mode

```
```
anvi-show-collections-and-bins -p PROFILE.db
```
OK. Let’s run the interactive interface again with the CONCOCT collection:
```
anvi-interactive -p PROFILE.db \
                 -c CONTIGS.db \
                 --collection-autoload CONCOCT

```
How do your bins compare?

Anvi’o also has a script called anvi-script-merge-collections to merge multiple files from external binning results into a single merged file (don’t ask why):
```
anvi-script-merge-collections -c CONTIGS.db \
                              -i additional-files/external-binning-results/*.txt \
                              -o collections.tsv
```
Visualize all binning results:
```
anvi-interactive -p PROFILE.db \
                 -c CONTIGS.db \
                 -A collections.tsv
```

To visually emphasize relationships between bins, the authors created a state file for us to match colors where bins match
```
anvi-import-state --state additional-files/state-files/state-merged.json \
                  --name default \
                  -p PROFILE.db
```
```
anvi-interactive -p PROFILE.db \
                 -c CONTIGS.db \
                 -A collections.tsv

```
Based on this, would you refine your bins? Probably, but for now let's keep the ones you refined. 


## 📝 Assignment due next class on Canvas
1. With std_coverage.txt file, make a plot of your choice to show chaning abundance of bins across the samples. Provide the code you use if using R. 
