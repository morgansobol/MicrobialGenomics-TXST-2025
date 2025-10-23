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

The last pre-processing step we need to do is remove low-counts and adjust for library size. We filter out genes with very low counts, because those genes:
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
