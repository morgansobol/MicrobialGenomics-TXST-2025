# Week XX: 16S Amplicon analysis - Part 2

In this tutorial, we will continue working with 16S amplicon data, this time comparing samples and performing statistical analysis. We will continue working in R, but handing off our Dada2 processed samples to the package Phyloseq. 

We will continue following Dada2's tutorial + Mike Lee's tutorial on Dada2, using Mike's data (with some modifications) below. Thanks all for the great documentation!
https://benjjneb.github.io/dada2/tutorial_1_8.html 
https://astrobiomike.github.io/amplicon/dada2_workflow_ex#differential-abundance-analysis-with-deseq2
https://blogs.oregonstate.edu/earthmotes/2021/09/28/dada2-pipeline-for-16s-datasets-in-r/

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:

- 

## 🧪 Step 1: Setting up the working environment and reading in the data
Open a new script to work in and start by loading dada2 and these other packages:
```R
library(dada2); packageVersion("dada2")
library(tidyverse) ; packageVersion("tidyverse") 
library(phyloseq) ; packageVersion("phyloseq") 
library(vegan) ; packageVersion("vegan") 
library(DESeq2) ; packageVersion("DESeq2")
library(dendextend) ; packageVersion("dendextend") 
library(viridis) ; packageVersion("viridis") 

setwd("[insert path to directory you want to be in]")

list.files() # make sure our files from last time are here
# ASVs-no-contam.fa
# ASVs_counts-no-contam.tsv
# ASVs_taxonomy-no-contam.tsv
# ASVs_counts-no-contam.tsv

# ok, now moving on

rm(list=ls()) # remove any prior objects, start from a clean slate

# Read in data ------------------------------------------------------------

# Here, blanks are removed, i.e. -c(1:4), because we already decontaminated the data
count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
             check.names=F, sep="\t")[ , -c(1:4)]

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
           row.names=1, check.names=F, sep="\t"))

```
Oops, we need a file Mike provided from the data_dir. Open the terminal in R and let's `cp` it here.
```bash
cp ../../data_dir/dada2/sample_info.tsv .
```

Ok, proceed with loading it:
```R
sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1,
                   check.names=F, sep="\t")
  
 # and setting the color column to be of type "character", which helps later
sample_info_tab$color <- as.character(sample_info_tab$color)

# take a peek at the data we have on each sample. 
# The rownames of the sample_info_tab must be the same as the column names of the count table for graphing later on. This should be the sample sites in the same order.

sample_info_tab
```

Taking a look at the `sample_info_tab`, we see it has the 16 samples as rows, and four columns: 
* 1) “temperature” for the temperature of the venting water was where collected;
  2) “type” indicating if that sample is a water sample, rock, or biofilm;
  3) a characteristics column called “char” that just serves to distinguish between the main types of rocks (glassy, altered, or carbonate);
  4) “color”, which has different R colors we can use later for plotting.

This table can be made anywhere (e.g. in R, in excel, at the command line), you just need to make sure you read it into R properly (which you should always check just like we did here to make sure it came in correctly).

## 🧪 Step 2: Make a Phyloseq Object
```R
# Make a Phyloseq object --------------------------------------------------

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
tax_tab_phy <- tax_table(tax_tab)
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
ASV_physeq
```

Should see something like this, telling us that all the info we want is in the Phyloseq object:
```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2498 taxa and 16 samples ]
sample_data() Sample Data:       [ 16 samples by 4 sample variables ]
tax_table()   Taxonomy Table:    [ 2498 taxa by 7 taxonomic ranks ]
```

## 🧪 Step 3: Normalize data

First, let's check the sequence abundance at all our sites and graph it for visualization. Sequence runs don't always go as planned, so we should expect some differences in sampling depths in our samples, which is most likely not truly reflective of that sample's true sequence abundance. 

```R
# Subsample data ----------------------------------------------------------

Sequences_per_sample <- sample_sums(phyloseq_object)
Sequences_per_sample
sample_richness <- data.frame(Sequences = Sequences_per_sample, Sample_Site = Group)
```

You see, we have almost a whole order of magnitude difference in sequence abundances. If we compare diversity metrics directly without accounting for these differences, the results could be biased toward the deeper sequenced samples, which may appear artificially more diverse simply because more sequences were recovered.

To make fair statistical comparisons across samples, we need to account for these differences in sequencing depth. Rarefaction does this by subsampling each sample down to the same number of sequences, repeated many times, and averaging across the subsamplings. This ensures that observed differences in diversity reflect biology rather than uneven sequencing effort

>[!NOTE]
>There has been some debate in the field about what normalization methods to apply.
>
> McMurdie & Holmes 2014 claim too much data is lost when you use rarefraction, and instead suggest to use Variance Stabilizing Transformation offered thorugh DEseq2. (paper here: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
>
>However, Dr. Pat Schloss at the University of Michichan, who is a microbial ecologist that wrote one of the OG 16S amplicon softwares (Mothur) and is very knowledgeable in this field, wrote a rebuttal to that paper and providing evidence that rarefraction/subsampling is still superior in most cases. (paper: https://journals.asm.org/doi/10.1128/msphere.00355-23?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed).
>
>He has also published a few Youtube videos on this and other stuff, I suggest checking it out if interested. Ex: https://www.youtube.com/watch?v=t5qXPIS-ECU&list=PLmNrK_nkqBpJuhS93PYC-Xr5oqur7IIWf&index=2&ab_channel=RiffomonasProject





DESeq2 creates negative values resulting from the log-like transformation are set to zero, causing the method to ignore many rare species completely. Additionally, many microbial environments are extremely variable in microbial composition, which would *violate* DESeq normalization assumptions of a constant abundance of a majority of species and of a balance of increased/decreased abundance for those species that do change.
```R
# rarefy data
asvs_rare = t(rrarefy(t(asvs_count), sample_size))
```
To get an idea of whether your rarefied data is representative of the samples, we can also plot a rarefaction curve:
```R
# plot curves
rarecurve(t(asvs_count), step=1000, label=F)
# add line to show rarefaction
abline(v=sample_size, col="red")
```
You can see the depth to which they have been rarefied by the red line. Ideally, the line intersects the sample curves when they are relatively horizontal.

## 🧪 Step 3: Calculate diversity
In ecology, the diversity of a community, or for us a sample, refers to various different measures. They can be classified as measures of either alpha diversity, which is calculated on a single sample, or beta diversity, which is calculated by comparing multiple samples.

### Beta diversity
Beta diversity involves calculating metrics such as distances or dissimilarities based on pairwise comparisons of samples so we can relate samples to each other. 

Typically the first thing to do is generate some exploratory visualizations like ordinations and hierarchical clusterings. These give you a quick overview of how your samples relate to each other and can be a way to check for problems like batch effects.

**Hierarchical clustering**
We’re going to use Euclidean distances {explain euclidiean distances) to generate some exploratory visualizations of our samples. Since differences in sampling depths between samples can influence distance/dissimilarity metrics, we first need to somehow normalize across our samples.
```R
euc_dist <- dist(t(asvs_rare))
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
plot(euc_clust) 
    # but i like to change them to dendrograms for two reasons:
        # 1) it's easier to color the dendrogram plot by groups
        # 2) if wanted you can rotate clusters with the rotate() 
        #    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")
```
So from our first peek, the broadest clusters separate the biofilm, carbonate, and water samples from the basalt rocks, which are the black and brown labels. And those form two distinct clusters, with samples R8–R11 separate from the others (R1-R6, and R12). R8-R11 (black) were all of the glassier type of basalt with thin (~1-2 mm), smooth exteriors, while the rest (R1-R6, and R12; brown) had more highly altered, thick (>1 cm) outer rinds (excluding the oddball carbonate which isn’t a basalt, R7). This is starting to suggest that level of alteration of the basalt may be correlated with community structure. If we look at the map figure again (below), we can also see that level of alteration also co-varies with whether samples were collected from the northern or southern end of the outcrop as all of the more highly altered basalts were collected from the northern end.

**Ordination**
```R
# making our phyloseq object with transformed table
rare_count_phy <- otu_table(asvs_rare, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
rare_physeq <- phyloseq(rare, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
rare_pcoa <- ordinate(rare_physeq, method="MDS", distance="euclidean")
eigen_vals <- rare_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(rare_physeq, rare_pcoa, color="char") + 
    geom_point(size=1) + labs(col="type") + 
    geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
    coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
    scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + 
    theme_bw() + theme(legend.position="none")
```
This is just providing us with a different overview of how our samples relate to each other. Inferences that are consistent with the hierarchical clustering above can be considered a bit more robust if the same general trends emerge from both approaches. It’s important to remember that these are exploratory visualizations and do not say anything statistically about our samples. But our initial exploration here shows us the rock microbial communities seem to be more similar to each other than they are to the water samples, and focusing on just the basalts (brown and black labels), these visualizations both suggest their communities may correlate with level of exterior alteration.

### Alpha diversity
Alpha diversity entails using summary metrics that describe individual samples, and it is a very tricky thing when working with amplicon data. If and when I use any alpha diversity metrics, I mostly consider them useful for relative comparisons of samples from the same experiment, and think of them as a just another summary metric (and not some sort of absolute truth about a sample).

**Rarefraction**

First thing’s first, it is *not* okay to use rarefaction curves to estimate total richness of a sample, or to extrapolate anything from them really, but they can still be useful in getting a sense of the landscape of your samples, depending on the data. We’ll be using the `rarecurve()` function from the package vegan here. Note that vegan expects rows to be samples and observations (our ASVs here) to be columns, which is why we transpose the first table in the command with `t()`.

```R
rarecurve(t(count_tab), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=F)

    # and adding a vertical line at the fewest seqs in any sample
abline(v=(min(rowSums(t(count_tab)))))
```

In this plot, samples are colored the same way as above, and the black vertical line represents the sampling depth of the sample with the least amount of sequences (a bottom water sample, BW1, in this case). This view suggests that the rock samples have a greater richness (unique number of sequences recovered) than the water samples or the biofilm sample – based on where they all cross the vertical line of lowest sampling depth, which is not necessarily predictive of where they’d end up had they been sampled to greater depth. And again, just focusing on the brown and black lines for the two types of basalts we have, they seem to show similar trends within their respective groups that suggest the more highly altered basalts (brown lines) may host more microbial communities with greater richness than the glassier basalts (black lines).

**Richness and diversity measurements

Two basic measures of a sample’s species diversity are:
* Richness - the number of species present in a sample
* Evenness - how equal the abundances of different species in a sample are**

Here we’re going to plot Chao1 richness esimates and Shannon diversity values. Chao1 is a richness estimator, “richness” being the total number of distinct units in our sample, “distinct units” being whatever we happen to be measuring (ASVs in our case here). And Shannon’s diversity index is a metric of diversity. The term diversity includes “richness” (the total number of our distinct units) and “evenness” (the relative proportions of all of our distinct units). Again, these are really just metrics to help contrast our samples within an experiment, and should not be considered “true” values of anything or be compared across studies.

We are going to go back to using the phyloseq package for this to use the function `plot_richness()`.
```R
# first we need to create a phyloseq object using our un-transformed count table
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

# and now we can call the plot_richness() function on our phyloseq object
plot_richness(ASV_physeq, color="char", measures=c("Chao1", "Shannon")) + 
    scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) +
    theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```
Before saying anything about this I’d like to stress again that these are not interpretable as “real” numbers of anything (due to the nature of amplicon data), but they can still be useful as relative metrics of comparison within a study. For example, we again see from this that the more highly altered basalts seem to host communities that are more diverse and have a higher richness than the other rocks, and that the water and biofilm samples are less diverse than the rocks.

And just for another quick example of why phyloseq is pretty awesome, let’s look at how easy it is to plot by grouping sample types while still coloring by the characteristics column:
```R
plot_richness(ASV_physeq, x="type", color="char", measures=c("Chao1", "Shannon")) + 
    scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) +
    theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```
One last note on interpretation here, don’t forget Chao1 is richness and Shannon is diversity, and what these mean as discussed above. Take for example the biofilm sample (green) above. It seems to have a higher estimated richness than the two water samples, but a lower Shannon diversity than both water samples. This suggests that the water samples likely have a greater “evenness”; or to put it another way, even though the biofilm may have more biological units (our ASVs here), it may be largely dominated by only a few of them.

## 🧪 Step 4: Taxonomic summaries
Don’t forget that the taxonomy called here was done rapidly and by default has to sacrifice some specificity for speed. For the sequences that become important in your story, you should absolutely pull them out and BLAST them, and possibly integrate them into phylogenetic trees to get a more robust idea of who they are most closely related to.

Here we’ll make some broad-level summarization figures. Phyloseq is also very useful for parsing things down by taxonomy now that we’ve got all that information in there. So I’ll be using that where I’m familiar with it, but unfortunately I’m not that familiar with probably about 90% of its functionality. So here you’ll see my ugly way of doing things, and then I’m counting on you to shoot me cleaner code so we can update this walkthrough for everyone’s benefit 🙂

Let’s make a summary of all major taxa proportions across all samples, then summaries for just the water samples and just the rock samples. To start, we need to parse our count matrix by taxonomy. How you want to break things down will depend on your data and your question, as usual. Here, we’ll just generate a table of proportions of each phylum, and then breakdown the Proteobacteria to the class level.

```R
# using phyloseq to make a count table that has summed all ASVs
        # that were in the same phylum
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="phylum")) 

  # making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="phylum"))[,"phylum"]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

    # we also have to account for sequences that weren't assigned any
    # taxonomy even at the phylum level 
    # these came into R as 'NAs' in the taxonomy table, but their counts are
    # still in the count table
    # so we can get that value for each sample by subtracting the column sums
    # of this new table (that has everything that had a phylum assigned to it)
    # from the column sums of the starting count table (that has all
    # representative sequences)
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
    # and we'll add this row to our phylum count table:
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

    # now we'll remove the Proteobacteria, so we can next add them back in
    # broken down by class
temp_major_taxa_counts_tab <- phyla_and_unidentified_counts_tab[!row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ]

    # making count table broken down by class (contains classes beyond the
    # Proteobacteria too at this point)
class_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="class")) 

    # making a table that holds the phylum and class level info
class_tax_phy_tab <- tax_table(tax_glom(ASV_physeq, taxrank="class")) 

phy_tmp_vec <- class_tax_phy_tab[,2]
class_tmp_vec <- class_tax_phy_tab[,3]
rows_tmp <- row.names(class_tax_phy_tab)
class_tax_tab <- data.frame("phylum"=phy_tmp_vec, "class"=class_tmp_vec, row.names = rows_tmp)

    # making a vector of just the Proteobacteria classes
proteo_classes_vec <- as.vector(class_tax_tab[class_tax_tab$phylum == "Proteobacteria", "class"])

    # changing the row names like above so that they correspond to the taxonomy,
    # rather than an ASV identifier
rownames(class_counts_tab) <- as.vector(class_tax_tab$class) 

    # making a table of the counts of the Proteobacterial classes
proteo_class_counts_tab <- class_counts_tab[row.names(class_counts_tab) %in% proteo_classes_vec, ] 

    # there are also possibly some some sequences that were resolved to the level
    # of Proteobacteria, but not any further, and therefore would be missing from
    # our class table
    # we can find the sum of them by subtracting the proteo class count table
    # from just the Proteobacteria row from the original phylum-level count table
proteo_no_class_annotated_counts <- phyla_and_unidentified_counts_tab[row.names(phyla_and_unidentified_counts_tab) %in% "Proteobacteria", ] - colSums(proteo_class_counts_tab)

  # now combining the tables:
major_taxa_counts_tab <- rbind(temp_major_taxa_counts_tab, proteo_class_counts_tab, "Unresolved_Proteobacteria"=proteo_no_class_annotated_counts)

    # and to check we didn't miss any other sequences, we can compare the column
    # sums to see if they are the same
    # if "TRUE", we know nothing fell through the cracks
identical(colSums(major_taxa_counts_tab), colSums(count_tab)) 

    # now we'll generate a proportions table for summarizing:
major_taxa_proportions_tab <- apply(major_taxa_counts_tab, 2, function(x) x/sum(x)*100)

  # if we check the dimensions of this table at this point
dim(major_taxa_proportions_tab)
    # we see there are currently 42 rows, which might be a little busy for a
    # summary figure
    # many of these taxa make up a very small percentage, so we're going to
    # filter some out
    # this is a completely arbitrary decision solely to ease visualization and
    # intepretation, entirely up to your data and you
    # here, we'll only keep rows (taxa) that make up greater than 5% in any
    # individual sample
temp_filt_major_taxa_proportions_tab <- data.frame(major_taxa_proportions_tab[apply(major_taxa_proportions_tab, 1, max) > 5, ])
  # checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
    # now we have 12, much more manageable for an overview figure

    # though each of the filtered taxa made up less than 5% alone, together they
    # may add up and should still be included in the overall summary
    # so we're going to add a row called "Other" that keeps track of how much we
    # filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(major_taxa_proportions_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

    ## don't worry if the numbers or taxonomy vary a little, this might happen due to different versions being used 
    ## from when this was initially put together
```
Now that we have a nice proportions table ready to go, we can make some figures with it. While not always all that informative, expecially at the level of resolution we’re using here (phyla and proteo classes only), we’ll make some stacked bar charts, boxplots, and some pie charts. We’ll use ggplot2 to do this, and for these types of plots it seems to be easiest to work with tables in narrow format. We’ll see what that means, how to transform the table, and then add some information for the samples to help with plotting.

```R
# first let's make a copy of our table that's safe for manipulating
filt_major_taxa_proportions_tab_for_plot <- filt_major_taxa_proportions_tab

    # and add a column of the taxa names so that it is within the table, rather
    # than just as row names (this makes working with ggplot easier)
filt_major_taxa_proportions_tab_for_plot$Major_Taxa <- row.names(filt_major_taxa_proportions_tab_for_plot)

    # now we'll transform the table into narrow, or long, format (also makes
    # plotting easier)
filt_major_taxa_proportions_tab_for_plot.g <- pivot_longer(filt_major_taxa_proportions_tab_for_plot, !Major_Taxa, names_to = "Sample", values_to = "Proportion") %>% data.frame()

    # take a look at the new table and compare it with the old one
head(filt_major_taxa_proportions_tab_for_plot.g)
head(filt_major_taxa_proportions_tab_for_plot)
    # manipulating tables like this is something you may need to do frequently in R

    # now we want a table with "color" and "characteristics" of each sample to
    # merge into our plotting table so we can use that more easily in our plotting
    # function
    # here we're making a new table by pulling what we want from the sample
    # information table
sample_info_for_merge<-data.frame("Sample"=row.names(sample_info_tab), "char"=sample_info_tab$char, "color"=sample_info_tab$color, stringsAsFactors=F)

    # and here we are merging this table with the plotting table we just made
    # (this is an awesome function!)
filt_major_taxa_proportions_tab_for_plot.g2 <- merge(filt_major_taxa_proportions_tab_for_plot.g, sample_info_for_merge)

    # and now we're ready to make some summary figures with our wonderfully
    # constructed table

    ## a good color scheme can be hard to find, i included the viridis package
    ## here because it's color-blind friendly and sometimes it's been really
    ## helpful for me, though this is not demonstrated in all of the following :/ 

    # one common way to look at this is with stacked bar charts for each taxon per sample:
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(x=Sample, y=Proportion, fill=Major_Taxa)) +
    geom_bar(width=0.6, stat="identity") +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, vjust=0.4, hjust=1), legend.title=element_blank()) +
    labs(x="Sample", y="% of 16S rRNA gene copies recovered", title="All samples")
```
Ok, that’s not helpful really at all in this case, but it’s here for the code example. No one likes stacked taxonomy barcharts, but we as humans just keep making them ¯\_(ツ)_/¯ Another way to look would be using boxplots where each box is a major taxon, with each point being colored based on its sample type.
```R
ggplot(filt_major_taxa_proportions_tab_for_plot.g2, aes(Major_Taxa, Proportion)) +
    geom_jitter(aes(color=factor(char), shape=factor(char)), size=2, width=0.15, height=0) +
    scale_color_manual(values=unique(filt_major_taxa_proportions_tab_for_plot.g2$color[order(filt_major_taxa_proportions_tab_for_plot.g2$char)])) +
    geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.title=element_blank()) +
    labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="All samples")
```
Don’t forget to keep in mind again that this was a very coarse level of resolution as we are using taxonomic classifications at the phylum and class ranks. This could partly be why things may look more similar between the rocks and water samples than we might expect, and why when looking at the ASV level – like we did with the exploratory visualizations above – we can see more clearly that these seem to host distinct communities. But let’s look at this for a second anyway. 

The biofilm sample (green triangles) clearly stands out as having the greatest proportion of seqs classified as coming from Alphaproteobacteria and those that were Unclassified. Three of the four “glassy” basalts (black plus signs) seem to have the greatest proportion of Gammaproteobacteria-derived sequences. And Cyanos, Desulfobacterota, and Firmicutes for the most part only seem to show up in one of the water samples. Another way to look at this would be to plot the water and rock samples separately, which might help tighten up some taxa boxplots if they have a different distribution between the two sample types.

    # let's set some helpful variables first:
bw_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$type == "water"]
rock_sample_IDs <- row.names(sample_info_tab)[sample_info_tab$type == "rock"]

    # first we need to subset our plotting table to include just the rock samples to plot
filt_major_taxa_proportions_rocks_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% rock_sample_IDs, ]
    # and then just the water samples
filt_major_taxa_proportions_water_samples_only_tab_for_plot.g <- filt_major_taxa_proportions_tab_for_plot.g2[filt_major_taxa_proportions_tab_for_plot.g2$Sample %in% bw_sample_IDs, ]
```R
    # and now we can use the same code as above just with whatever minor alterations we want
    # rock samples
ggplot(filt_major_taxa_proportions_rocks_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
    scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
    geom_jitter(aes(color=factor(char), shape=factor(char)), size=2, width=0.15, height=0) +
    scale_color_manual(values=unique(filt_major_taxa_proportions_rocks_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_rocks_only_tab_for_plot.g$char)])) +
    geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="top", legend.title=element_blank()) + # moved legend to top 
    labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Rock samples only")

    # water samples
ggplot(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g, aes(Major_Taxa, Proportion)) +
    scale_y_continuous(limits=c(0,50)) + # adding a setting for the y axis range so the rock and water plots are on the same scale
    geom_jitter(aes(color=factor(char)), size=2, width=0.15, height=0) +
    scale_color_manual(values=unique(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g$color[order(filt_major_taxa_proportions_water_samples_only_tab_for_plot.g$char)])) +
    geom_boxplot(fill=NA, outlier.color=NA) + theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
    labs(x="Major Taxa", y="% of 16S rRNA gene copies recovered", title="Bottom-water samples only")
    
    # probably don't need the boxplots for the 2 water samples, but let's be crazy
```
This shows us more clearly for instance that one of the two water samples had ~15% of its recovered 16S sequences classified as Firmicutes, while none of the 13 rock samples had more than 1%. It’s good of course to break down our data and look at it in all the different ways you can, this was just to demonstrate one example.

Last taxonomic summary we’ll go through in will just be some pie charts. This is mostly because I think it’s worth showing an example of using the `pivot_wider()` function to return our tables into “wide” format.

```R
    # notice we're leaving off the "char" and "color" columns, in the code, and be sure to peak at the tables after making them
rock_sample_major_taxa_proportion_tab <- filt_major_taxa_proportions_rocks_only_tab_for_plot.g[, c(1:3)] %>% pivot_wider(names_from = Major_Taxa, values_from = Proportion) %>% column_to_rownames("Sample") %>% t() %>% data.frame()
water_sample_major_taxa_proportion_tab <- filt_major_taxa_proportions_water_samples_only_tab_for_plot.g[, c(1:3)] %>% pivot_wider(names_from = Major_Taxa, values_from = Proportion) %>% column_to_rownames("Sample") %>% t() %>% data.frame()

    # summing each taxa across all samples for both groups 
rock_sample_summed_major_taxa_proportions_vec <- rowSums(rock_sample_major_taxa_proportion_tab)
water_sample_summed_major_taxa_proportions_vec <- rowSums(water_sample_major_taxa_proportion_tab)

rock_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(rock_sample_summed_major_taxa_proportions_vec), "Proportion"=rock_sample_summed_major_taxa_proportions_vec, row.names=NULL)
water_sample_major_taxa_summary_tab <- data.frame("Major_Taxa"=names(water_sample_summed_major_taxa_proportions_vec), "Proportion"=water_sample_summed_major_taxa_proportions_vec, row.names=NULL)

    # plotting just rocks
ggplot(data.frame(rock_sample_major_taxa_summary_tab), aes(x="Rock samples", y=Proportion, fill=Major_Taxa)) + 
    geom_bar(width=1, stat="identity") +
    coord_polar("y") +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Rock samples only") +
    theme_void() +
    theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())

    # and plotting just water samples
ggplot(data.frame(water_sample_major_taxa_summary_tab), aes(x="Bottom water samples", y=Proportion, fill=Major_Taxa)) + 
    geom_bar(width=1, stat="identity") +
    coord_polar("y") +
    scale_fill_viridis(discrete=TRUE) +
    ggtitle("Water samples only") +
    theme_void() +
    theme(plot.title = element_text(hjust=0.5), legend.title=element_blank())
```
Again, not very useful here. But here is how you might parse your dataset down by taxonomy to whatever level actually is useful and make some standard visualizations.

## 🧪 Step 5: Differential abundances
**Betadisper and permutational ANOVA**
As we saw earlier, we have some information about our samples in our sample info table. There are many ways to incorporate this information, but one of the first I typically go to is a permutational ANOVA test to see if any of the available information is indicative of community structure. Here we are going to test if there is a statistically signficant difference between our sample types. One way to do this is with the `betadisper` and `adonis` functions from the vegan package. `adonis` can tell us if there is a statistical difference between groups, but it has an assumption that must be met that we first need to check with `betadisper`, and that is that there is a sufficient level of homogeneity {explain} of dispersion within groups. If there is not, then `adonis` can be unreliable.
```R
anova(betadisper(euc_dist, sample_info_tab$type)) # 0.002
    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together, and due to variation in the permutations
```
Checking by all sample types, we get a significant result (0.002) from the betadisper test. This tells us that there is a difference between group dispersions, which means that we can’t trust the results of an adonis (permutational anova) test on this, because the assumption of homogenous within-group disperions is not met. This isn’t all that surprising considering how different the water and biofilm samples are from the rocks.

But a more interesting and specific question is “Do the rocks differ based on their level of exterior alteration?” So let’s try this looking at just the basalt rocks, based on their characteristics of glassy and altered.
```R
 # first we'll need to go back to our transformed table, and generate a
    # distance matrix only incorporating the basalt samples
    # and to help with that I'm making a variable that holds all basalt rock
    # names (just removing the single calcium carbonate sample, R7)
basalt_sample_IDs <- rock_sample_IDs[!rock_sample_IDs %in% "R7"]

    # new distance matrix of only basalts
basalt_euc_dist <- dist(t(vst_trans_count_tab[ , colnames(vst_trans_count_tab) %in% basalt_sample_IDs]))

    # and now making a sample info table with just the basalts
basalt_sample_info_tab <- sample_info_tab[row.names(sample_info_tab) %in% basalt_sample_IDs, ]

    # running betadisper on just these based on level of alteration as shown in the images above:
anova(betadisper(basalt_euc_dist, basalt_sample_info_tab$char)) # 0.7
    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together, and due to variation in the permutations
```
Looking at just the two basalt groups, glassy vs the more highly altered with thick outer rinds, we do not find a significant difference between their within-group dispersions (0.7). So we can now test if the groups host statistically different communities based on adonis, having met this assumption.

```R
adonis(basalt_euc_dist~basalt_sample_info_tab$char) # 0.003
    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together, and due to variation in the permutations
```
And with a significance level of 0.003, this gives us our first statistical evidence that there is actually a difference in microbial communities hosted by the more highly altered basalts as compared to the glassier less altered basalts, pretty cool!

It can be useful to incorporate this statistic into a visualization. So let’s make a new PCoA of just the basalts, and slap our proud significance on there.
```R
   # making our phyloseq object with transformed table
basalt_vst_count_phy <- otu_table(vst_trans_count_tab[, colnames(vst_trans_count_tab) %in% basalt_sample_IDs], taxa_are_rows=T)
basalt_sample_info_tab_phy <- sample_data(basalt_sample_info_tab)
basalt_vst_physeq <- phyloseq(basalt_vst_count_phy, basalt_sample_info_tab_phy)

    # generating and visualizing the PCoA with phyloseq
basalt_vst_pcoa <- ordinate(basalt_vst_physeq, method="MDS", distance="euclidean")
basalt_eigen_vals <- basalt_vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

    # and making our new ordination of just basalts with our adonis statistic
plot_ordination(basalt_vst_physeq, basalt_vst_pcoa, color="char") + 
    labs(col="type") + geom_point(size=1) + 
    geom_text(aes(label=rownames(basalt_sample_info_tab), hjust=0.3, vjust=-0.4)) + 
    annotate("text", x=25, y=68, label="Highly altered vs glassy") +
    annotate("text", x=25, y=62, label="Permutational ANOVA = 0.003") + 
    coord_fixed(sqrt(basalt_eigen_vals[2]/basalt_eigen_vals[1])) + ggtitle("PCoA - basalts only") + 
    scale_color_manual(values=unique(basalt_sample_info_tab$color[order(basalt_sample_info_tab$char)])) + 
    theme_bw() + theme(legend.position="none")
```
**Differential abundance analysis with DESeq2**
First, it’s important to keep in mind that:

Recovered 16S rRNA gene copy numbers do not equal organism abundance.

That said, recovered 16S rRNA gene copy numbers do represent… well, numbers of recovered 16S rRNA gene copies. So long as you’re interpreting them that way, and thinking of your system in the appropriate way, you can perform differential abundance testing to test for which representative sequences have significantly different copy-number counts between samples – which can be useful information and guide the generation of hypotheses. One tool that can be used for this is DESeq2, which we used above to transform our count table for beta diversity plots.

Now that we’ve found a statistical difference between our two rock groups, this is one way we can try to find out which ASVs (and possibly which taxa) are contributing to that difference. If you are going to use DESeq2, be sure to carefully go over their thorough manual and other information you can find here, particularly this very informative page.

We are going to take advantage of another phyloseq convenience, and use the phyloseq_to_deseq2 function to make our DESeq2 object.
```R
   # first making a basalt-only phyloseq object of non-transformed values (as that is what DESeq2 operates on
basalt_count_phy <- otu_table(count_tab[, colnames(count_tab) %in% basalt_sample_IDs], taxa_are_rows=T)
basalt_count_physeq <- phyloseq(basalt_count_phy, basalt_sample_info_tab_phy)
  
    # now converting our phyloseq object to a deseq object
basalt_deseq <- phyloseq_to_deseq2(basalt_count_physeq, ~char)

    # and running deseq standard analysis:
basalt_deseq <- DESeq(basalt_deseq)
```
The DESeq() function is doing a lot of things. Be sure you look into it and understand conceptually what is going on here, it is well detailed in their manual.

We can now access the results. In our setup here, we only have 2 groups, so what is being contrasted is pretty straightforward. Normally, you will tell the results() function which groups you would like to be contrasted (all were done at the DESeq2() function call, but we would parse those initial results by specifying when we pull them out with the results() function). We will also provide the p-value we wish to use to filter the results with later, as recommended by the ?results help page, with the “alpha” argument.

```R
 # pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_altered_vs_glassy <- results(basalt_deseq, alpha=0.01, contrast=c("char", "altered", "glassy"))

    # we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_altered_vs_glassy) 
    # this tells us out of ~1,800 ASVs, with adj-p < 0.01, there are 7 increased when comparing altered basalts to glassy basalts, and about 6 decreased
    # "decreased" in this case means at a lower count abundance in the altered basalts than in the glassy basalts, and "increased" means greater proportion in altered than in glassy
    # remember, this is done with a drastically reduced dataset, which is hindering the capabilities here quite a bit i'm sure

    # let's subset this table to only include these that pass our specified significance level
sigtab_res_deseq_altered_vs_glassy <- deseq_res_altered_vs_glassy[which(deseq_res_altered_vs_glassy$padj < 0.01), ]

    # now we can see this table only contains those we consider significantly differentially abundant
summary(sigtab_res_deseq_altered_vs_glassy) 

    # next let's stitch that together with these ASV's taxonomic annotations for a quick look at both together
sigtab_deseq_altered_vs_glassy_with_tax <- cbind(as(sigtab_res_deseq_altered_vs_glassy, "data.frame"), as(tax_table(ASV_physeq)[row.names(sigtab_res_deseq_altered_vs_glassy), ], "matrix"))

    # and now let's sort that table by the baseMean column
sigtab_deseq_altered_vs_glassy_with_tax[order(sigtab_deseq_altered_vs_glassy_with_tax$baseMean, decreasing=T), ]

    # this puts a sequence derived from a Rhizobiales at the second to highest (first is unclassified) that was detected in ~7 log2fold greater abundance in the glassy basalts than in the highly altered basalts
```
If you glance through the taxonomy of our significant table here, you’ll see several have the same designations. It’s possible this is one of those cases where the single-nucleotide resolution approach more inhibits our cause than helps it. You can imagine that with organisms having multiple copies of the 16S rRNA gene, which may not be identical, this could be muddying what we’re looking for here by splitting the signal up and weaking it. Another way to look at this would be to sum the ASVs by the same genus designations, or to go back and cluster them into some form of OTU (after identifying ASVs) – in which case we’d still be using the ASV units, but then clustering them at some arbitrary level to see if that level of resolution is more revealing for the system we’re looking at.

## Bonus: Environmental drivers of community diversity

## 📝 Assignment due next class on Canvas
Perform a top-level exploration of one of the sample outputs and answer these three questions. 
1. What are the dominant taxa in the sample you chose?
2. Tell me something about what this taxa is known for.
3. What % of your sample is unclassified? (I.e. not assignment to the genus level). 

