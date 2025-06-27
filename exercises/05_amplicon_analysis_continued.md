# Week XX: 16S Amplicon analysis - Part 2

In this tutorial we will continue working with 16S amplicon data, this time comparing samples and performing statistical analysis. We will continue working in R. 

We will continue following Dada2's tutorial + Mike Lee's tutorial on Dada2, using Mike's data (with some modifications) below. Thanks Mike!

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

setwd("~/dada2_amplicon_ex_workflow")

list.files() # make sure our files from last time are here
# ASVs-no-contam.fa
# ASVs_counts-no-contam.tsv
# ASVs_taxonomy-no-contam.tsv
# ASVs_counts-no-contam.tsv

# ok, now moving on

rm(list=ls())
  
count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
             check.names=F, sep="\t")[ , -c(1:4)]   ## remove blanks from count table

```
>[!Note]
>Since we’ve already used decontam to remove likely contaminants, we’re dropping the “blank” samples from our count table which are the first 4 columns. That’s what’s being done by the [ , -c(1:4)] part at the end there.

```R
tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
           row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.table("sample_info.tsv", header=T, row.names=1,
                   check.names=F, sep="\t")
  
    # and setting the color column to be of type "character", which helps later
sample_info_tab$color <- as.character(sample_info_tab$color)

sample_info_tab # to take a peek at the data we have on each sample
```
Taking a look at the `sample_info_tab`, we see it has the 16 samples as rows, and four columns: 
* 1) “temperature” for the temperature of the venting water was where collected;
  2) “type” indicating if that sample is a blank, water sample, rock, or biofilm;
  3) a characteristics column called “char” that just serves to distinguish between the main types of rocks (glassy,       altered, or carbonate);
  4) 4) “color”, which has different R colors we can use later for plotting.

This table can be made anywhere (e.g. in R, in excel, at the command line), you just need to make sure you read it into R properly (which you should always check just like we did here to make sure it came in correctly).

## 🧪 Step 2: Normalize data
Why do we need to normalize our data? Well sequence runs don't always go as planned, so we should expect some differences in sampling depths in our samples, that may not be reflective of that sample's true sequence abundance. 

There has been some debate in the field about what normalization method to apply. Ultimately, it comes down to your sample types. Here we will briefly dicuss the main methods, rarefraction and 

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
You can see the depth to which you have rarefied by the red line. Ideally the line is intersecting the sample curves when they are relatively horizontal.

## 🧪 Step 3: Calculate diversity
In ecology, the diversity of a community, or for us a sample, refers to various different measures. They can be classified as measures of either alpha diversity, which is calculated on a single sample, or beta diversity, which is calculated by comparing multiple samples.



