# Week 04: 16S Amplicon analysis - Part 2

In this tutorial, we will continue working with 16S amplicon data, this time comparing samples and performing statistical analysis. We will continue working in R, but handing off our Dada2 processed samples to the package Phyloseq. 

We will continue following Dada2's tutorial + Mike Lee's tutorial on Dada2, using Mike's data (with some modifications) below. Thanks all for the great documentation!
https://benjjneb.github.io/dada2/tutorial_1_8.html 
https://astrobiomike.github.io/amplicon/dada2_workflow_ex#differential-abundance-analysis-with-deseq2
https://blogs.oregonstate.edu/earthmotes/2021/09/28/dada2-pipeline-for-16s-datasets-in-r/

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:

- Understand why normalization is important for amplicon data.
- Visualize and discuss community differences using diversity metrics.
- Interpret diversity indices to describe within and between sample differences.

## 🧪 Step 1: Setting up the working environment and reading in the data
Open a new script to work in and start by loading dada2 and these other packages:
```R
install.packages("remotes")
remotes::install_github("cpauvert/psadd")
library("psadd")

# If you are not already in the dada2 directory, setwd to get there
setwd("[insert path to directory you want to be in]")

list.files() # make sure our files from last time are here

# ok, now moving on


# Read in data ------------------------------------------------------------

# blanks are removed, i.e. -c(1:4), because we already decontaminated the data
# this also creates ASVs with no counts in any other sample, so we need to remove them
count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1, check.names=F, sep="\t")[ , -c(1:4)]
count_tab <- count_tab[rowSums(count_tab) > 0, ]

tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T, row.names=1, check.names=F, sep="\t"))

```

We will also need a file Mike provided from the data_dir. Open the terminal in R and let's `cp` it here.
```bash
cp ../../data_dir/sample_info.tsv .
```

Ok, proceed with loading it:
```R
sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1, check.names=F, sep="\t")
  
# Setting the color column to be of type "character", which helps later
sample_info_tab$color <- as.character(sample_info_tab$color)

# take a peek at the data we have on each sample.

sample_info_tab

# The rownames of the sample_info_tab must be the same as the column names of the count table for graphing later on. This should be the sample sites in the same order.

count_tab

```

Taking a look at the `sample_info_tab`, we see it has the 16 samples as rows, and four columns: 
  1) “temperature” for the temperature of the venting water was where collected;
  2) “type” indicating if that sample is a water sample, rock, or biofilm;
  3) a characteristics column called “char” that just serves to distinguish between the main types of rocks (glassy, altered, or carbonate);
  4) “color”, which has different R colors we can use later for plotting.

This table can be made anywhere (e.g. in R, in excel, at the command line), you just need to make sure you read it into R properly (which you should always check just like we did here to make sure it came in correctly).

## 🧪 Step 2: Make a Phyloseq Object
Phyloseq is an R package to import, store, analyze, and graphically display complex phylogenetic sequencing data that has already been clustered into OTUs or ASVs. It leverages and builds upon many of the tools available in R for ecology and phylogenetic analysis (vegan, ade, ape), while also using advanced/flexible graphic systems (ggplot2) to easily produce publication-quality graphics. Check out more, including tutorials on what else it can do: https://joey711.github.io/phyloseq/ 

It stores sequencing data as a single, self-consistent, self-describing experiment-level object, making it easier to use. Let's make that "object". 

```R
# Make a Phyloseq object --------------------------------------------------

count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
tax_tab_phy <- tax_table(tax_tab)
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
ASV_physeq
```

You should see something like this, telling us that all the info we want is in the Phyloseq object:
```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2498 taxa and 16 samples ]
sample_data() Sample Data:       [ 16 samples by 4 sample variables ]
tax_table()   Taxonomy Table:    [ 2498 taxa by 7 taxonomic ranks ]
```

## 🧪 Step 3: Normalize data

First, let's check the sequence abundance at all our sites and graph it for visualization. Sequence runs don't always go as planned, so we should expect some differences in sampling depths in our samples, which is most likely not truly reflective of that sample's true sequence abundance. 

```R
# Normalize data ----------------------------------------------------------
seqs_per_sample <- sample_sums(ASV_physeq)
seqs_per_sample

# make that into a dataframe
sample_richness <- data.frame(Sequences = seqs_per_sample, Sample_Site = "Group")
```

You see, we have almost a whole order of magnitude difference in sequence abundances. If we compare diversity metrics directly without accounting for these differences, the results could be biased toward the deeper sequenced samples, which may appear artificially more diverse simply because more sequences were recovered.

To make fair statistical comparisons across samples, we need to account for these differences in sequencing depth. Rarefaction does this by subsampling each sample down to the same number of sequences, repeated many times, and averaging across the subsamplings. This ensures that observed differences in diversity reflect biology rather than uneven sequencing effort

> [!NOTE]
> There has been some debate in the field about what normalization methods to apply.
>
> McMurdie & Holmes 2014 claim that too much data is lost when you use rarefaction, and instead suggest using the Variance Stabilizing Transformation offered through DESeq2. But many microbial environments are extremely variable in microbial composition, which would *violate* DESeq normalization assumptions of a constant abundance of a majority of species and of a balance of increased/decreased abundance for those species that do change. (paper here: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531)
>
>Indeed, Dr. Pat Schloss at the University of Michigan, a microbial ecologist who wrote one of the OG 16S amplicon softwares (Mothur) and is very knowledgeable in this field, wrote a rebuttal to that paper, providing evidence that *true* rarefaction/subsampling is still superior in dealing with uneven sequence depth. (paper 1: https://journals.asm.org/doi/10.1128/msphere.00355-23?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed, paper 2: https://journals.asm.org/doi/full/10.1128/msphere.00354-23) 
>
>He has also published a few Youtube videos on this and other stuff, I suggest checking it out if interested. https://www.youtube.com/watch?v=t5qXPIS-ECU&list=PLmNrK_nkqBpJuhS93PYC-Xr5oqur7IIWf&index=2&ab_channel=RiffomonasProject

Before using rarefaction (with multiple sampling!) to normalize our data, let's look at a rarefaction curve to get a sense of how sequencing depth relates to observed diversity of ASVs.

```R
# rarefraction curve 
rarecurve(t(count_tab), step=100, col=sample_info_tab$color, lwd=2, ylab="ASVs", label=T, legend=TRUE)
abline(v=(min(rowSums(t(count_tab)))))

```
This view suggests that the rock samples have a greater richness (unique number of sequences recovered) than the water samples or the biofilm sample – based on where they all cross the vertical line of lowest sampling depth, which is not necessarily predictive of where they’d end up had they been sampled to greater depth. And again, just focusing on the brown and black lines for the two types of basalts we have, they seem to show similar trends within their respective groups that suggest the more highly altered basalts (brown lines) may host more microbial communities with greater richness than the glassier basalts (black lines).

Note the line we drew with `abline` shows us what would happen if we just subsampled all samples only once at the minimum values of sequences in a sample, i.e. 1897. We would lose some diversity that way, to the point of McMurdie & Holmes. Let's test this out. We will do rarefaction with and without multiple sampling.

```R
# set normalization depth as the lowest number of sequences in a sample
depth_target = min(sample_richness)
depth_target

# first, without multiple subsampling
# rngseed = 123, fixes the starting point of the random number generator, ensuring that the subsampling of 1897 sequences will always return the same subset of sequences each time the code is run.
asv_rarefy <- rarefy_even_depth(ASV_physeq, sample.size=depth_target, rngseed=123)

# says we lost 356 ASVs. Let's confirm that another way.
ntaxa(ASV_physeq) # original dataset before normalization
ntaxa(asv_rarefy) 

```
So, yes, 356 ASVs were lost. Let's see what happens when we average across multiple subsamplings.
We will need to install the package first called `psadd`, which is an extension to phyloseq. 

```R

# same target depth, but now we want it to repeat the subsampling 1000 times. We also set a "seed" of 123 (random number), which fixes the starting point of the random number generator, ensuring that the subsampling of 1897 sequences will always return the same subset of sequences each time the code is run.

asv_multi_rare <- multiple_rrarefy(ASV_physeq, sample.size=depth_target, 1000, 123)

# confirm it worked first
seqs_per_sample <- sample_sums(asv_multi_rare)
seqs_per_sample

# Let's check if we lost any taxa. 
ntaxa(ASV_physeq)
ntaxa(asv_multi_rare)

# write that to a new file
asv_norm <- as(otu_table(asv_multi_rare), "matrix")
write.table(asv_norm, "ASVs_normalized.tsv", sep="\t", quote=F, col.names=NA)

```

We did not lose any taxa with this type of normalization, so I feel more confident about this approach. 

One con to this normalization is that subsampling creates non-integer numbers (i.e., not whole numbers) and this is an issue for some diversity indices (which we will calculate later) that require integers.
We need to round ASVs with counts between 0-1 to 1 instead of 0, otherwise we will lose the rare taxa in the sample. Then we will round the other numbers as is. Is this the right thing to do? Depends on your data.
For your sake, I tested dropping counts <0. The diversity results did not significantly change, what did happen is that we would have lost >100 ASVs total. So I feel its right not round <0 counts to 1. All we can do is clearly state what we did to be transparent. 🤷‍♀️

Let's do that rounding using a function in R:
```R
# Extract ASV matrix
asv_mat <- as(otu_table(asv_multi_rare), "matrix")

# Custom rounding: bump (0,1) up to 1, round everything else
asv_mat_round <- apply(asv_mat, c(1,2), function(x) {
  if (x > 0 & x < 1) {
    return(1L)
  } else {
    return(round(x,0))
  }
})

asv_norm_round <- data.frame(otu_table(asv_mat_round, taxa_are_rows = taxa_are_rows(asv_multi_rare)))

write.table(asv_norm_round, "ASVs_counts_rounded.tsv",
            sep="\t", quote=F, col.names=NA)

# Let's see how that changed our sequence counts
sample_sums(asv_norm_round)
min(sample_sums(asv_norm_round))
max(sample_sums(asv_norm_round))
# we "inflate" some samples with rare counts

# Rebuild phyloseq object with new count table
ASV_physeq <- phyloseq(
  asv_norm_round,
  tax_table(tax_tab_phy),
  sample_data(sample_info_tab_phy)
)
```

Now we can continue our analysis. 

## 🧪 Step 4: Plotting ASV abundance

You've probably seen one of those barplot figures in a paper somewhere that shows how the richness and evenness of taxa differ between samples.
Let's make one first before we assess any stats. We will need to turn our ASV counts into relative abundances, i.e., dividing each ASV’s count in a sample by the total number of sequences in that sample so that the values represent proportions that sum to one.

```R
# Abundance Analysis ------------------------------------------------------

asv_rel <- transform_sample_counts(ASV_physeq, function(x) x / sum(x))

# check it worked. Should see that every sample sum is equal to 1
colSums(otu_table(asv_rel)) %>% head()

# prepare data for plotting

# 1. Turn the phyloseq object into a numeric matrix for R to use
asv_mat_rel <- as(otu_table(asv_rel), "matrix")

# 2. Collapse ASVs to the phylum level
asv_phylum <- tax_glom(asv_rel, "Phylum", NArm=FALSE)
View(asv_phylum)

# 3. “Melt” the phyloseq object into a long-format data frame for ggplot.
df <- psmelt(asv_phylum)
View(df)

# 4. Make the plot
library(ggplot2)

ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# default colors are awful, let's change that

library(RColorBrewer)

n <- length(unique(df$Phylum))   # number of phyla in your data
mycols <- colorRampPalette(brewer.pal(12, "Paired"))(n)

ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Now, for better visualization, it could be beneficial to collapse rare taxa (<1% abundant) into their own category "Other", like so:
threshold <- 0.01  # collapse anything <1% abundance 
df$Phylum2 <- df$Phylum
df$Phylum2[df$Abundance < threshold] <- "Other"

n <- length(unique(df$Phylum2))
mycols <- colorRampPalette(brewer.pal(12, "Paired"))(n)

ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum2)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = mycols) +
  theme_bw()

```
## 🧪 Step 5: Calculate alpha diversity

A diversity index is a quantitative measure that is used to assess the level of diversity or variety within a particular system, such as a microbial community. 

**Alpha diversity** describes the diversity within a single community/sample. It considers the number of different species in that sample (also referred to as species richness). Additionally, it can take the abundance of each species into account to measure how evenly taxa are distributed across the sample (also referred to as species evenness). 

* Chao1 = Measures total _richness_, so observed + rare taxa inferred from singletons/doubletons. The higher the value, the greater the number of ASVs. 
* Shannon's = Measures _diversity_ as both richness and evenness (relative proportions of our ASVs). The value increases as you add more taxa, even if they are rare.
* Simpson's = also measures both richness and evenness, but weights more on dominant taxa and is less sensitive to rare ones. 1 = very even; close to 0 = one or a few taxa dominate.

To compare alpha diversity across samples, would be to ask if the mean or median of these calculated indices differs across groups.

> [WARNING!]
> These are just some metrics to help compare & contrast our samples within an experiment, and should **not** be considered “true” values of any ASV. 

```R
# Alpha diversity ---------------------------------------------------------

# We call on phyloseq's estimate_richness() function to calculate alpha diversity using Chao1, Shannon, and Simpson
asv_alpha <- estimate_richness(ASV_physeq, measures = c("Observed", "Chao1", "Shannon", "Simpson"))

# add sample info and save output
asv_alpha <- cbind(asv_alpha, as(sample_data(ASV_physeq), "data.frame"))
write.table(asv_alpha, "ASVs_alpha-diversity.tsv", sep="\t", quote=F, col.names=NA)

# We call on phyloseq's plot_richness() function on our phyloseq object to plot the data
plot_richness(ASV_physeq , color="char", measures=c("Observed", "Chao1", "Shannon", "Simpson")) + scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + theme_bw() + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
```
What can we say about alpha diversity across the samples?

We can also ask, are samples significantly different based on sample type or sample alteration ("char")? Microbiome data is compositional and generally violates many assumptions in statistical analyses, especially that of normality, so nonparametric tests (i.e., those that do not rely on assumptions about the distribution of the data) are often used. Moreover, given the number of pairwise comparisons that must be made to compare the samples, a multiple-testing correction should be applied to the P-value.

Kruskal–Wallis and Wilcoxon rank-sum are rank-based statistical tests. They require each group to have at least two observations so ranks can be meaningfully compared. 
> Wilcoxon compares two groups, i.e.  _Is diversity different between rock vs water samples?_
> Pairwise-wilcoxon compares groups to one another, i.e. _Does diversity differ between water vs glassy, water vs altered, or glassy vs altered?_
> Kruskal compares multiple groups, i.e. _Is diversity different across water, glassy, and altered rocks?_

```R
# Now let's subset the data to remove biolfim and carbonate since they are only one sample

asv_sub <- subset(asv_alpha, char %in% c("water","glassy","altered"))
asv_sub

# Test differences Chao1
kruskal.test(Chao1 ~ char, data=asv_sub)

pairwise.wilcox.test(
  asv_sub$Chao1,
  asv_sub$char,
  p.adjust.method = "BH")
```
Where do we see significant differences? 
What about for Shannon and Simpson?

## 🧪 Step 6: Calculate beta diversity

**Beta diversity**, also called "between-sample diversity", is a measurement of the distance, or difference, between samples. It involves calculating metrics such as distances or dissimilarities based on pairwise comparisons of samples so we can relate samples to each other. What stats do for us is determine the overall variation in the distance matrix and test whether groups of samples differ in community composition. 

Typically, you would generate some exploratory visualizations like ordinations and hierarchical clusterings for an overview of how your samples relate to each other. 
Let's look at all the different distance approaches Phyloseq can use:
```R
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
```
We’re going to use Bray-Curtis dissimilarity to cluster samples that are similar to one another based on ASV profiles. 
Bray-curtis looks at shared abundance between two samples. 
> If samples share many taxa with similar abundances → low dissimilarity (close to 0).
> If they have very different taxa or very different abundances → high dissimilarity (close to 1)

```R
asv_dist <- phyloseq::distance(ASV_physeq, method = "bray")
View(as.matrix(asv_dist))
```
**Hierarchical clustering**

```R
# dendogram
# method = "average", means the distance between two clusters is defined as the average of all pairwise distances between the samples in cluster X and cluster Y.
asv_hclust <- hclust(asv_dist, method="average")
asv_dend <- as.dendrogram(asv_hclust, hang=0.1)

# plot dendogram and color by char type
# Make sure rownames(sample_info_tab) are your sample IDs
color_map <- setNames(as.character(sample_info_tab$color),
                      rownames(sample_info_tab))

# Reorder colors according to the dendrogram labels
dend_cols <- color_map[labels(asv_dend)]

# Apply colors
labels_colors(asv_dend) <- dend_cols

plot(asv_dend, ylab="Bray-Curtis Distance")
```
So from our first peek, the broadest clusters separate the biofilm, carbonate, glassy, and water samples from the altered basalt rocks, which are the brown labels. R8-R11 (black) were all of the glassier type of basalt with thin (~1-2 mm), smooth exteriors, while the rest (R1-R6, and R12; brown) had more highly altered, thick (>1 cm) outer rinds (excluding the oddball carbonate which isn’t a basalt, R7). 

**Ordination clustering**
```R
asv_pcoa <- ordinate(ASV_physeq, method="PCoA", distance="bray")

plot_ordination(ASV_physeq, asv_pcoa, color="char") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(sample_info_tab$color[order(sample_info_tab$char)])) + theme_bw() + theme(legend.position="none")
```

**Further statistical testing**

The permutational multivariate analysis of variance (PERMANOVA) test is frequently used in comparing beta diversity, i.e. whether community composition differs among groups. It is an ANOVA for distance matrices. Although nonparametric, it assumes homogeneity of variability. If variability differs strongly between groups, significant results could be misleading. 

One way to do this is with the `betadisper` and `adonis` functions from the `vegan` package. `adonis` can tell us if there is a statistical difference between groups, but it has an assumption that must be met that we first need to check with `betadisper`, and that is that there is a sufficient level of homogeneity of variability (dispersion) within groups. If there is not, then `adonis` can be unreliable.

```R
anova(betadisper(asv_dist, sample_info_tab$type))
# So at least one type (water vs rock vs biofilm) is more variable than the others.

adonis2(asv_dist~sample_info_tab$type)

anova(betadisper(asv_dist, sample_info_tab$char))
# Some “char” groups are tightly clustered, others very spread out from one another.

adonis2(asv_dist~sample_info_tab$char)

```

Since `adonis` and `betadisper` are both significant, interpret cautiously: some groups may simply be more variable. That itself can be biologically meaningful. If one environment (say, biofilm) has communities that are very inconsistent, that’s an interesting ecological result, not necessarily an “error.” I think the "issue" is _the_ biofilm sample since it is so different from the others. It is not that the data is wrong, just a statistical concern for comparison. 
We would just report the results from `betadisper` and `adonis`


## Finale 

Our data altogether is starting to suggest that the level of alteration of the basalt may be correlated with community structure. If we look at the map figure again (below), we can also see that level of alteration also co-varies with whether samples were collected from the northern or southern end of the outcrop as all of the more highly altered basalts were collected from the northern end.

<img width="800" height="436" alt="image" src="https://github.com/user-attachments/assets/f3aed91d-1f8e-4d12-827f-ed3b5edb9a55" />


> [!TIP]
> After all that is said and done, you should know that Phyloseq offers _Shiny-Phyloseq_, which is a web-browser GUI to where you can point and click instead of write code to do your analysis and make figures. Of course, this does not replace the flexibility of coding your own data in R, but could be useful as a first exploration of your data.
>To try it out, you need to first install in R:
```R
install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")
```
## 

## 📝 Assignment due next class on Canvas
Perform a top-level exploration of the sample outputs and answer this question. 
1. Determine the _prevalence_ of ASVs across all samples. Here we will define as the number of samples in which an ASV appears at least once, i.e. how many samples each ASV is found in. Was there a single ASV that had the highest "prevalence", or were there multiple at equal prevalence? If so, what were they/it? 



