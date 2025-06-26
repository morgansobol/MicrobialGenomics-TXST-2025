# Week XX: 16S Amplicon analysis - Part 1

Welcome to Week XX! 
This week, you will finally start analyzing data, specially performing 16S amplicon analysis. From here on out, unless noted, we are working in R, not at the Unix-like command line. If you find you need some refresher in R basics check out this page: https://astrobiomike.github.io/R/ 

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:

- 

## Quick refresher on commonly used R scripts 

Where am I?
```R
getwd()
```
Change my working directory to the desired location
```R
setwd("[insert path to directory you want to be in]")
```
 List the contents of the directory
```R
dir()
```
Import a csv file
```R
iris_df <-read.csv("path/to/a/dataset/called/iris.csv",
                   header = TRUE, sep=",", strip.white=TRUE, stringsAsFactors=FALSE)
```
List column and row names of the data
```R
colnames(iris_df)
rownames(iris_df)
```
See the structure of the data
```R
str(iris_df)
```
View the data 
```R
view(iris_df)
```
Get basics stats
```R
summary(iris_df)
```

## 🧪 Exercise 1: Setting up the working environment
```R
library(dada2)
packageVersion("dada2") # 1.11.5 when this was initially put together, though might be different in the binder or conda installation, that's ok!

setwd("~/dada2_amplicon_ex_workflow")

list.files() # make sure what we think is here is actually here

## first we're setting a few variables we're going to use ##
  # one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("samples", what="character")

  # one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_sub_R1_trimmed.fq.gz")
  # and one with the reverse
reverse_reads <- paste0(samples, "_sub_R2_trimmed.fq.gz")

  # and variables holding file names for the forward and reverse
  # filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")
```
## 🧪 Exercise 2: More quality trimming/filtering
We did a filtering step above with cutadapt (where we eliminated reads that had imperfect or missing primers and those that were shorter than 215 bps or longer than 285), but in DADA2 we’ll implement a trimming step as well (where we trim reads down based on some quality threshold rather than throwing the read away). 

Since we’re potentially shortening reads further, we’re again going to include another minimum-length filtering component. We can also take advantage of a handy quality plotting function that DADA2 provides to visualize how you’re reads are doing, plotQualityProfile(). By running that on our variables that hold all of our forward and reverse read filenames, we can easily generate plots for all samples or for a subset of them. So let’s take a peak at that to help decide our trimming lengths:

```R
plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)
  # and just plotting the last 4 samples of the reverse reads
plotQualityProfile(reverse_reads[17:20])
```

All forwards look pretty similar to eachother, and all reverses look pretty similar to eachother, but worse than the forwards, which is common – chemistry gets tired 😞

Here’s the output of the last four samples’ reverse reads:
[insert pic]

On these plots, the bases are along the x-axis, and the quality score on the y-axis. The black underlying heatmap shows the frequency of each score at each base position, the green line is the mean quality score at that base position, the orange is the median, and the dashed orange lines show the quartiles. The red line at the bottom shows what percentage of reads are that length

>[!NOTE]
>In Phred talk the difference between a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100. In this case, since we have full overlap with these primers and the sequencing performed (515f-806r, 2x300), we can be pretty conservative and trim these up a bit more. But it’s important to think about your primers and the overlap you’re going to have.
>
>Here, our primers span 515-806 (291 bases), and we cut off the primers which were 39 bps total, so we are expecting to span, nominally, 252 bases. If we trimmed forward and reverse here down to 100 bps each, we would not span those 252 bases and this would cause problems later because we won’t be able to merge our forward and reverse reads.
>
>Make sure you’re considering this based on your data. Here, I’m going to cut the forward reads at 250 and the reverse reads at 200 – roughly where both sets maintain a median quality of 30 or above – and then see how things look. But we also want to set a minimum length to filter out those that are too short to overlap (by default, this function truncates reads at the first instance of a quality score of 2, this is how we could end up with reads shorter than what we are explicitly trimming them down to).

In DADA2, this quality-filtering step is done with the filterAndTrim() function:
```R
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                rm.phix=TRUE, minLen=175, truncLen=c(250,200))
```
Here, the first and third arguments (“forward_reads” and “reverse_reads”) are the variables holding our input files, which are our primer-trimmed output fastq files from cutadapt. The second and fourth are the variables holding the file names of the output forward and reverse seqs from this function. And then we have a few parameters explicitly specified:
- maxEE is the quality filtering threshold being applied based on the expected errors and in this case we are saying we want to throw the read away if it is likely to have more than 2 erroneous base calls (we are specifying for both the forward and reverse reads separately).
- rm.phix removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring. And minLen is setting the minimum length reads we want to keep after trimming.
- As mentioned above, the trimming occurring beyond what we set with truncLen is coming from a default setting, truncQ, which is set to 2 unless we specify otherwise, meaning it trims all bases after the first quality score of 2 it comes across in a read.
- There is also an additional filtering default parameter that is removing any sequences containing any Ns, maxN, set to 0 by default.
- Then we have our truncLen parameter setting the minimum size to trim the forward and reverse reads to in order to keep the quality scores roughly above 30 overall.

As mentioned, the output read files were named in those variables we made above (“filtered_forward_reads” and “filtered_reverse_reads”), so those files were created when we ran the function – we can see them if we run list.files() in R, or by checking in our working directory in the terminal:

Let's check the object we generated called filtered_out, containing how many reads went in and how many reads made it out. How do we view this object?

[insert pic]

Now plot the reads like before using the plotQualityProfile function, *but* make sure you do it for the filtered reads this time. How do we do that? 

[insert pic]

## 🧪 Exercise 3: Generating an error model of our data
The dada2 tool that we will use for this analysis uses a statistic approach that aims to predict if sequences are actually true biological sequences or are partly generated by errors that appear during sequencing.

Each sequencing run, even when all goes well, will have its own subtle variations to its error profile. This step tries to assess that for both the forward and reverse reads. It is one of the more computationally intensive steps of the workflow. We will use the multithread=TRUE

