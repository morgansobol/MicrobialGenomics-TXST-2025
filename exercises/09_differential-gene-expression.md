## Week 9: Differential Expression Analysis

Studies of differential gene expression generally aim to determine the differences between a baseline condition and an experimental condition. In the simplest case, expression levels are compared between two conditions—for example, control vs. treatment or mutant vs. wild‑type. More complex experiments may include additional experimental factors, such as multiple mutants or different doses of a drug.

In this tutorial, you will compare RNA‑Seq data in order to detect genes that are expressed with statistically significant differences. 

We will be working in R this time.

## 🧠 Learning Objectives

By the end of this exercise, you should be able to:
* Gain more familiarity with R. 
* Understand the importance of normalization of data to make samples comparable.
* Quantify gene expression differences between experimental conditions.
* Identify significantly differentially expressed genes. 


## 🧪 Step 1: Reading in the data and install initial software
Like the last few weeks, you need to make a new `week9` directory inside your `microgenomics-2025` directory. 
Make a `work_dir` inside `week9`. Then, from Canvas Module 9, download the data_dir.zip file and use `unzip` to unzip the directory. 
Print your `week9` path using `pwd`. 

Now open up R Studio and start a new script. To set your path in R:
```R
setwd("copy/path/from/pwd")
```
Use the command `getwd()` to verify that `setwd()` worked. 

We will first use the `edgeR` package for DGE analysis. 
Publication here: https://academic.oup.com/nar/article/53/2/gkaf018/7973897 

`BiocManager` is a tool to install software for biological analysis in R. Install it first in case it's not already installed. Then we will use it to install `edgeR`. 
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.21")

# load BiocManager
library(BiocManager)

BiocManager::install("edgeR")

# load edgeR
library(edgeR)

# load ggplot2 for plotting
library(ggplot2)

# load reshape to help us reformat data matrices
library(reshape2)

```

## 🧪 Step 2: Exploration and pre-processing of our data

For this tutorial, we will consider samples from two different conditions: a microbe grown with fructose vs. with glycerol as a carbon source. How do they behave differently? 

The count matrix was created by mapping sequence reads against a reference genome. Rows are genes/transcripts, and columns are samples. For each condition, we have three replicate samples. To read the data into R:

```R
countTable <- read.table("data_dir/matrix.count.txt",
                         header = TRUE, sep = "	", row.names = 1)
head(countTable)

```
RNA-seq counts vary widely across genes within samples. To reduce this variation and aid visualization, we apply a log2 transformation, i.e. compress large counts and expand small ones. 

However, you can't take the log of 0, so we add a count of +1 before the log2 transformation to ensure genes with zero counts are included and we don't end up with zero count data. 

```r
pseudoCount <- log2(countTable + 1)
```

Let's plot our non-transformed and transformed counts to see what I mean. 
Our data is currently in "wide" format (samples as columns), but ggplot2 uses "long" format (samples as rows). 
We will use reshapes 'melt' function to convert our matrix. 

```R
count_melt = melt(countTable, variable.name = "Samples",id.vars=NULL)
pseudoCount_melt = melt(pseudoCount, variable.name = "Samples",id.vars=NULL)
```

Additionlly, we want to create a "condition" column for plotting later that informs which treatment samples belong to. 

* The regular expression "_.$" means “underscore + one last character (the replicate number)”.
* sub("_.$", "", x) removes that part, leaving just "Fructose" or "Glycerol".
```R
count_melt = data.frame(count_melt, Condition = sub("_.$","",count_melt$Samples))
pseudoCount_melt = data.frame(pseudoCount_melt, Condition = sub("_.$","",pseudoCount_melt$Samples))
```
This creates a new column "Condition" that groups replicates correctly.

Now we can plot our data using boxplots to show the variation in counts across samples. 
This quick QC step tells you whether:

* Library sizes or sequencing depths differ a lot between samples.
* There are potential outliers (e.g., one replicate has a shifted distribution).

```R
ggplot(count_melt, aes(x = Samples, y = value, fill = Condition)) + 
  geom_boxplot() + ylab("Counts") + theme_light()

ggplot(pseudoCount_melt, aes(x = Samples, y = value, fill = Condition)) + 
  geom_boxplot() + ylab("Counts") + theme_light()
```
As you can see the wide variation in counts makes it hard to visualize. Let's stick with the `pseuduoCount" data for visuals. 

All boxes look roughly aligned (similar medians and spread), so normalization will likely perform well, and your DGE results will be trustworthy (hopefully). 

We can also plot samples on a 2D plot to see how they cluster based on similarity of gene counts. We want to ensure that the treatments and their replicates cluster together. 
```R
# generate the dissimilarity matrix
dataMDS <- plotMDS(pseudoCount, plot = FALSE)

# Use sample names from your matrix
sample <- colnames(pseudoCount)

# Build a data frame with dissimilarity matrix and additional data
dfMDS <- data.frame(
  x = dataMDS$x,
  y = dataMDS$y,
  sample = sample,
  condition = sub("_.$", "", sample)
)

# plot
ggplot(dfMDS, aes(x, y, color = condition)) +
  geom_point(show.legend = FALSE) +
  theme_bw() +
  labs(title = "MDS Plot", x = "Dimension 1", y = "Dimension 2") +
  theme(plot.title = element_text(size = 10, hjust = 0.5)) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3, show.legend = FALSE)
```

The last pre-processing step we need to do is to remove low-counts and adjust for library size. We filter out genes with very low counts, because those genes:
* are often not truly expressed (just background noise),
* and they make the statistical model unstable (too many zeros).

It’s not appropriate to apply a filter on top of the pseudoCount matrix if we are asking whether a gene is measurably expressed. 

We will use "counts per million", or "cpm" to scale raw counts to reads per million total reads. In other words, out of every million reads in this sample, how many came from gene X? This lets us compare expression levels fairly between samples, while preserving magnitude differences. 
```R
cpmsTable <- cpm(countTable)
head(cpmsTable)
```

Now filter out genes with CPM ≤ 1 in fewer than "n" samples, where "n" is the size of the smallest group of replicates, i.e. "Keep genes that have CPM > 1 in at least 3 samples".
```R
keep <- rowSums(cpmsTable > 1) >= 3
countTableFilter <- countTable[keep, ]
removed <- countTable[!keep, ]

dim(countTableFilter)  # kept
dim(removed)           # removed
dim(countTable)        # original
```

Ok, stop here and we will continue DGE later in lab. 

## 🧪 Step 3: DGE analysis

We will create a `DGEList` to store our filtered count information and sample/condition info.
```R
# make sure names are as intended
names(countTableFilter)

# create the condition vector derived from the names
conditions <- sub("_.$", "", names(countTableFilter))

# create the dgelist
dge <- DGEList(counts = countTableFilter, group = conditions)
dge
# Note: norm.factors == 1 indicates not yet normalized
```

Now we are ready to normalize our data for differential analysis. `edgeR` uses the command `calcNormFactors`. Let's take a look at what our options are for this command using the HELP tab. 

The recommended normalization method in `edgeR` is TMM, or Trimmed Mean of M values. This takes into account sequencing depth and sample count composition. TMM normalization adjusts library sizes based on the assumption that most genes are not differentially expressed. It corrects for library composition bias (not just library size), e.g. if one sample has a few genes with extremely high counts, without correction, many other genes might appear artificially low in that library. The trimming + mean calculation helps mitigate that. And we don't need to account for gene length in this case since we are comparing the same gene across different samples. 

Let's break down what is going on:

1. Pick a reference sample to normalize all samples against.
For all samples, the 75th quartile is calculated on based on counts per million scaling, which separates the top 25% of genes with the highest counts from the lower 75%. The chosen reference sample is the one whose 75th quartile value is closest to the average 75th quartile of scaled counts across _all_ libraries.

<img width="1154" height="339" alt="image" src="https://github.com/user-attachments/assets/0296394c-cfc7-4186-a98d-167440b4afaa" />

2. Compare genes from each sample to the reference to create a scaling factor.
The scaling factor in TMM normalization is a number that tells `edgeR` how much to rescale that sample so that the majority of counts for the genes line up with the reference sample. Extreme high or low counts get trimmed. This makes expression values comparable across samples.

4. Apply the scaling factor to adjust library sizes.
Each sample’s original library size is multiplied by its TMM scaling factor to create an effective library size.
These effective sizes are then used to compute normalized counts per million (CPM).
This step corrects for differences in sequencing depth and sample composition, ensuring that any observed differences in gene expression reflect true biological variation, not technical bias.

<img width="1000" height="500" alt="image" src="https://github.com/user-attachments/assets/35d18593-fa2d-4262-9648-14818e989294" />

If this is difficult to digest, I suggest watching this Youtube video: https://www.youtube.com/watch?v=Wdt6jdi-NQo 

Now all of that happens in this one command:

```R
dge <- calcNormFactors(dge)
dge

# Re‑check MDS after normalization to see if samples still cluster by condition
plotMDS(dge, col = c(2,2,2,3,3,3)) 
```

Ok, to test for differential expression, `edgeR` needs to know how much variation is expected for each gene.

```R
# estimates biological variability for the genes across replicates
dge <- estimateDisp(dge)

# Here we ask for a single global estimate of variability across all genes ("baseline")
dge$common.dispersion   

```

Biological Coefficient of Variation (BCV) = √0.0168 ≈ 0.13 (13%)

That means, on average, gene expression levels vary by about 13% between biological replicates after accounting for sequencing depth and normalization.

**Now** we can test the hypothesis, i.e. if genes are expressed differently between conditions. 

The `exactTest()` function asks, for every gene, “Are the counts in Fructose and Glycerol different enough that it’s unlikely to be just random variation?”

It uses:
* Your normalized counts (after calcNormFactors)
* The dispersion estimates (from estimateDisp)
* Your group labels (stored in dge$samples$group)

to test the null hypothesis:

𝐻0: mean expression of a gene is the same in both groups

vs. 

𝐻1: a gene is differentially expressed between groups.

By default, `edgeR` uses the first factor level as the baseline.
```R
dge$samples$group
[1] Fructose Fructose Fructose Glycerol Glycerol Glycerol
Levels: Fructose Glycerol
```
So, we will be asking if genes in Glycerol are different from Fructose treatment.  

Run the exact test:
```R
et <- exactTest(dge)
```
This give you an output with the following info:
* logFC
* logCPM
* PValue
* FDR (False Discovery Rate; Benjamini–Hochberg corrected p-value)

Let's summarize the DGEs.
```R
 # top 10 by default
topTags(et)
 # all genes   
resEdgeR <- topTags(et, n = Inf)$table
```

We will select DGEs to analyze based on some thresholds typically used in DGE analysis:
* FDR ≤ 0.01 = You expect fewer than 1% false positives among DEGs.
* |logFC| ≥ 1.5 = DEGs should have a fold change greater than +1.5 or lower than -1.5

```R
genesDE <- decideTests(et, lfc = 1.5, p.value = 0.01)
summary(genesDE)
```

We can also just directly filter out table using those thresholds
```R
gDE <- resEdgeR[abs(resEdgeR$logFC) >= 1.5 & resEdgeR$FDR <= 0.01, ]
dim(gDE)
```

Now let's label those results with meaningful names:
```R
resEdgeR[rownames(genesDE), "Expression"] <- data.frame(genesDE[, "Glycerol-Fructose"])

resEdgeR[resEdgeR$Expression ==  1, "Expression"] <- "UpInGlycerol"
resEdgeR[resEdgeR$Expression == -1, "Expression"] <- "UpInFructose"
resEdgeR[resEdgeR$Expression ==  0, "Expression"] <- "non-DE"
```

Now save results:
```R
write.table(resEdgeR, "work_dir/All_Results_EdgeR.txt", sep = "	", quote = FALSE)
write.table(genesDE, "work_dir/Genes_DE_edgeR_calls.txt", sep = "	", quote = FALSE)
```

## 🧪 Step 4: Plot results

Great, now we can plot our results, starting with the volcano plot. We use −log10 (FDR) in a volcano plot purely for visualization, not for any statistical calculation, because FDR values are very small. 

```R
resEdgeR$negLogFDR <- -log10(resEdgeR$FDR)

ggplot(resEdgeR, aes(x = logFC, y = negLogFDR, color = Expression)) +
  geom_point(shape = 1) + theme_bw() +
  scale_color_manual(breaks = c("UpInGlycerol", "non-DE", "UpInFructose"),
                     values = c("red", "black", "green")) +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "blue", size = 0.75) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue", size = 0.75)

ggsave("work_dir/volcano_plot.pdf")
```

On the other hand, we can plot a heatmap:
```R
# install ComplexHeatmap
BiocManager::install("ComplexHeatmap")

# load it
library(ComplexHeatmap)

# Normalized counts
countNormalized <- cpm(dge, normalized.lib.size = TRUE)

# Select genes classified as DE
DEGenes <- rownames(resEdgeR[resEdgeR$Expression != "non-DE", ])
signifGenNorm <- countNormalized[DEGenes, ]

# Z‑score across genes (rows)
data4heatmap <- t(scale(t(signifGenNorm), center = TRUE, scale = TRUE))

Heatmap(
  data4heatmap,
  name = "z-score",
  col = topo.colors(12),
  show_row_names = TRUE
)
```



