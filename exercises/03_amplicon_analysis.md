# Week 2: 16S Amplicon analysis - Part 1

Welcome to Week 2! 
This week, you will finally start analyzing data. Specifically, we will be performing 16S amplicon analysis. From here on out, unless noted, we are working in R, not at the Unix-like command line. If you find you need a refresher in R basics, check out these pages: 
 * https://astrobiomike.github.io/R/.
 * http://r-tutorial.nl/
 * https://rstudio-education.github.io/hopr/
 * Google is your friend too...


We will be following DADA2's tutorial + Mike Lee's tutorial on DADA2, using Mike's data (with some modifications) below. Thanks Mike!
https://benjjneb.github.io/dada2/tutorial.html

https://astrobiomike.github.io/amplicon/dada2_workflow_ex

---
## 🧠 Learning Objectives

By the end of this exercise, you should be able to:

- Interpret DADA2 error learning and ASV inference, and describe why exact sequence variants (ASVs) are preferred over OTUs.
- Assess read merging and chimera removal, including identifying potential causes of sequence loss at these stages.
- Classify ASVs taxonomically and recognize limitations of reference databases. 
- Learn the importance of sampling blanks and decontamination

## Commonly used R scripts 

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
filename_df <-read.csv("path/to/a/dataset/ending/in/filename.csv",
                   header = TRUE, sep=",", strip.white=TRUE, stringsAsFactors=FALSE)
```
List column and row names of the data
```R
colnames(filename_df)
rownames(filename_df)
```
See the structure of the data
```R
str(filename_df)
```
View the data 
```R
view(filename_df)
```
Get basic stats
```R
summary(filename_df)
```

Write a table/data matrix as a tab-delimited file
```R
write.table(object, "filename.tsv", sep="\t")
```

Write a new variable for R to use. A variable is a name that stores a value or some data.  
```bash
new_df = filename_df
new_df <- filename_df
filename_df -> new_df
```
All three mean the same thing: you are creating a new variable called `new_df` that stores whatever is in `filename_df` (for example, a dataframe).

Ok, let's begin!

First, let's learn a bit of where the sequences we have been working on come from:

Mike and his team were exploring an underwater mountain ~3 km down at the bottom of the Pacific Ocean that serves as a low-temperature (~5-10°C) hydrothermal venting site. This amplicon dataset was generated from DNA extracted from crushed basalts collected from across the mountain with the goal being to begin characterizing the microbial communities of these deep-sea rocks. No one had ever been here before, so as is often the purpose of marker-gene sequencing, this was just a broad-level community survey. The sequencing was done on the Illumina MiSeq platform with 2x300 paired-end sequencing using primers targeting the V4 region of the 16S rRNA gene. There are 20 samples total: 4 extraction “blanks” (nothing added to DNA extraction kit), 2 bottom-water samples, 13 rocks, and one biofilm scraped off a rock. 

In the following figure, overlain on the map are the rock sample collection locations, and the panes on the right show examples of the 3 distinct types of rocks collected: 1) basalts with highly altered, thick outer rinds (>1 cm); 2) basalts that were smooth, glassy, thin exteriors (~1-2 mm); and 3) one calcified carbonate.

<img width="800" height="436" alt="image" src="https://github.com/user-attachments/assets/aa197211-0392-4b63-998a-6e32de69efb5" />

This work was published, so I encourage you to check it out: https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2015.01470/full 

## 🧪 Step 1: Setting up the working environment
```R
##### Script for analysis with DADA2 ####

# Installing Packages -----------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

# During installation, it will ask:
# Do you want to install from sources the package which needs compilation? (yes/no/cancel) type yes, press enter
# Update all/some/none? Type "a" press enter

#restart R (not sure why we have to do this, but oh well)
.rs.restartR()

install.packages("Rcpp")
install.packages(c("httpuv","later","promises"))
install.packages("tidyverse")
install.packages("dendextend")
install.packages("viridis")

# Again, it will ask, Update all/some/none? Type "a" press enter

# Now run each line separately:
# If it asks, # Do you want to install from sources the package which needs compilation? (yes/no/cancel) type yes, press enter
# Or, Update all/some/none?, Type "a", press enter

BiocManager::install("decontam")
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
BiocManager::install("vegan")

# load the packages
library(dada2); packageVersion("dada2") #1.34.0
library(decontam); packageVersion("decontam") #1.26.0
library(DECIPHER); packageVersion("decipher") #3.2.0
library(tidyverse) ; packageVersion("tidyverse") # 1.3.1
library(phyloseq) ; packageVersion("phyloseq") # 1.50.0
library(vegan) ; packageVersion("vegan") # 2.7.1
library(dendextend) ; packageVersion("dendextend") # 1.19.1
library(viridis) ; packageVersion("viridis") # 0.6.5

# Setting Our Path --------------------------------------------------------

getwd() # where you are currently

# On Mac, you can navigate to the folder where your trimmed read files are, select the 3 dots as the top, select Get Info, copy the info following "Where", then past below. 

setwd("[insert path to directory you want to be in]")

list.files() # make sure what we think is here is actually here
dir() # this works too

rm(list=ls()) # remove any prior objects, start from a clean slate

# Set Up Our Variables ----------------------------------------------------

## first we're setting a few variables we're going to use ##
  # one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("samples.txt", what="character")

  # one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_sub_R1_trimmed.fq.gz")
  # and one with the reverse
reverse_reads <- paste0(samples, "_sub_R2_trimmed.fq.gz")

  # and variables holding file names for the forward and reverse
  # filtered reads we're going to generate below
filtered_forward_reads <- paste0(samples, "_sub_R1_filtered.fq.gz")
filtered_reverse_reads <- paste0(samples, "_sub_R2_filtered.fq.gz")
```

## 🧪 Step 2: More quality trimming/filtering
We did a filtering step above with cutadapt (where we eliminated reads that had imperfect or missing primers and those that were shorter than 215 bps or longer than 285), but in DADA2 we’ll implement a trimming step as well (where we trim reads down based on some quality threshold rather than throwing the read away). 

Since we’re potentially shortening reads further, we’re again going to include another minimum-length filtering component. We can also take advantage of a handy quality plotting function that DADA2 provides to visualize how your reads are doing, `plotQualityProfile()`. By running that on our variables that hold all of our forward and reverse read filenames, we can easily generate plots for all samples or for a subset of them. So let’s take a peak at that to help decide our trimming lengths:

```R
# More QC -----------------------------------------------------------------

plotQualityProfile(forward_reads)
plotQualityProfile(reverse_reads)
  # and just plotting the last 4 samples of the reverse reads to get a closer look
plotQualityProfile(reverse_reads[17:20])
```

On these plots, the bases are along the x-axis, and the quality score on the y-axis. The black underlying heatmap shows the frequency of each score at each base position, the green line is the mean quality score at that base position, the orange is the median, and the dashed orange lines show the quartiles. The red line at the bottom shows what percentage of reads are that length

All forwards look pretty similar to eachother, and all reverses look pretty similar to eachother, but worse than the forwards, which is common – chemistry gets tired 😞

>[!NOTE]
>In Phred talk the difference between a quality score of 40 and a quality score of 20 is an expected error rate of 1 in 10,000 vs 1 in 100. In this case, since we have full overlap with these primers and the sequencing performed (515f-806r, 2x300), we can be pretty conservative and trim these up a bit more. But it’s important to think about your primers and the overlap you’re going to have.
>
>Here, our primers span 515-806 (291 bases), and we cut off the primers which were 39 bps total, so we are expecting to span, nominally, 252 bases. If we trimmed forward and reverse here down to 100 bps each, we would not span those 252 bases and this would cause problems later because we won’t be able to merge our forward and reverse reads.
>
>Make sure you’re considering this based on your data. Here, I’m going to cut the forward reads at 250 and the reverse reads at 200 – roughly where both sets maintain a median quality of 30 or above – and then see how things look. But we also want to set a minimum length to filter out those that are too short to overlap (by default, this function truncates reads at the first instance of a quality score of 2, this is how we could end up with reads shorter than what we are explicitly trimming them down to).

>[!TIP]
> Zymo developed a tool called Figaro (I have not used it myself yet) that can help you choose DADA2 trimming parameters: https://github.com/Zymo-Research/figaro#figaro
> If time permits, perhaps we can try later, but for now we will continue with trimming.

In DADA2, this quality-filtering step is done with the `filterAndTrim()` function:
```R
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                rm.phix=TRUE, minLen=175, truncLen=c(250,200))
```
Here, the first and third arguments (“forward_reads” and “reverse_reads”) are the variables holding our input files, which are our primer-trimmed output fastq files from cutadapt. 

The second and fourth are the variables holding the file names of the output forward and reverse seqs from this function. 

And then we have a few parameters explicitly specified:
- `maxEE` is the quality filtering threshold being applied based on the expected errors and in this case we are saying we want to throw the read away if it is likely to have more than 2 erroneous base calls (we are specifying for both the forward and reverse reads separately).
- `rm.phix` removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring. And minLen is setting the minimum length reads we want to keep after trimming.
- As mentioned above, the trimming occurring beyond what we set with truncLen is coming from a default setting, `truncQ`, which is set to 2 unless we specify otherwise, meaning it trims all bases after the first quality score of 2 it comes across in a read.
- There is also an additional filtering default parameter that is removing any sequences containing any Ns, `maxN`= 0 by default.
- Then we have our `truncLen` parameter setting the minimum size to trim the forward and reverse reads to in order to keep the quality scores roughly above 30 overall.

As mentioned, the output read files were named in those variables we made above (“filtered_forward_reads” and “filtered_reverse_reads”), so those files were created when we ran the function – we can see them if we run `list.files()` in R, or by checking in our working directory in the terminal. 

Let's check the object we generated called `filtered_out`, containing how many reads went in and how many reads made it out. 

Now plot the reads like before using the `plotQualityProfile` function, *but* make sure you do it for the filtered reads this time. Try this yourself by typing the commands in the script and running them. They should look much better now. 

## 🧪 Step 3: Generating an error model of our data
The dada2 tool that we will use for this analysis uses a statistic approach that aims to predict if sequences are actually true biological sequences or are partly generated by errors that appear during sequencing.

Each sequencing run, even when all goes well, will have its own subtle variations to its error profile. This step tries to assess that for both the forward and reverse reads. It is one of the more computationally intensive steps of the workflow. We will use the `multithread=TRUE` to speed this up. 
```R
# Generate an Error Model -------------------------------------------------

err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
```
Now we plot these results to see what errors we have (ignore the warnings). 
```R
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
```
The red line is what is expected based on the quality score, the black line represents the estimate, and the black dots represent the observed. Generally speaking, you want the observed (black dots) to track well with the estimated (black line). So things look good and we can move on!

## 🧪 Step 4: Dereplication

Dereplication is a common step in many amplicon processing workflows. Instead of keeping 100 identical sequences and doing all downstream processing to all 100, you can keep/process one of them, and just attach the number 100 to it. When DADA2 dereplicates sequences, it also generates a new quality-score profile of each unique sequence based on the average quality scores of each base of all of the sequences that were replicates of it. 

The dereplication step is technically no longer listed as part of the standard dada2 tutorial, as it is performed by the dada() step if given filenames.
But, it can be lighter on memory requirements to run them as separate steps like done here, so it is left this way here for the sake of keeping things going. 

```R
# Dereplication -----------------------------------------------------------

derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples
```

## 🧪 Step 5: Inferring ASVs
We are now ready to apply the core algorithim, `dada`, and see what it was made to do, that is to do its best to infer true biological sequences. 💪

The dada2 tool will inspect every sequence and decide, based on the error model, if a sequence is a real biological sequence with no errors or a sequence that contains errors.

This step can be run on individual samples, which is the least computationally intensive manner, or on all samples together, which increases the function’s ability to resolve low-abundance ASVs. 
>Imagine Sample A has 10,000 copies of sequence Z, and Sample B has 1 copy of sequence Z. Sequence Z would likely be filtered out of Sample B even though it was a “true” singleton among perhaps thousands of spurious singletons we needed to remove.

Because running all samples together on large datasets can become impractical computationally, the developers also added a way to try to combine the best of both worlds they refer to as *pseudo-pooling*, which is explained very nicely here: https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling.

This basically provides a way to tell Sample B that sequence Z is legit. But it’s noted at the end of the pseudo-pooling page that this is not always the best way to go, and it may depend on your experimental design which is likely more appropriate for your data – as usual. There are no one-size-fits-all solutions in bioinformatics! 
```R
dada_forward <- dada(derep_forward, err=err_forward_reads, pool="pseudo", multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, pool="pseudo", multithread=TRUE)
```

## 🧪 Step 6: Merging forward and reverse reads
Now DADA2 merges the forward and reverse ASVs to reconstruct our full target amplicon requiring the overlapping region to be identical between the two. By default it requires that at least 12 bps overlap, but in our case the overlap should be much greater. If you remember above we trimmed the forward reads to 250 and the reverse to 200, and our primers were 515f–806r. After cutting off the primers we’re expecting a typical amplicon size of around 260 bases, so our typical overlap should be up around 190. However, that’s estimated based on E. coli 16S rRNA gene positions and very back-of-the-envelope-esque of course, so to allow for true biological variation and such I’m going ot set the minimum overlap for this dataset for 170. I’m also setting the trimOverhang option to TRUE in case any of our reads go passed their opposite primers (which I wouldn’t expect based on our trimming, but is possible due to the region and sequencing method).

```R
# Merging Reads -----------------------------------------------------------

merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                    derep_reverse, trimOverhang=TRUE, minOverlap=170)

    # this object holds a lot of information that may be the first place you'd want to look if you want to start poking under the hood
class(merged_amplicons) # list
length(merged_amplicons) # 20 elements in this list, one for each of our samples
names(merged_amplicons) # the names() function gives us the name of each element of the list 

class(merged_amplicons$B1) # each element of the list is a dataframe that can be accessed and manipulated like any ordinary dataframe

names(merged_amplicons$B1) # the names() function on a dataframe gives you the column names
# "sequence"  "abundance" "forward"   "reverse"   "nmatch"    "nmismatch" "nindel"    "prefer"    "accept"
```

## 🧪 Step 7: Generating a count table
Now we can generate a count table with the makeSequenceTable() function. This is one of the main outputs from processing an amplicon dataset. You may have also heard this referred to as a biome table, or an OTU matrix.
```R
# Generate a Count Table --------------------------------------------------

seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2521

    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together
```

The count table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. 

## 🧪 Step 8: Remove chimeras
Chimeras are separate, individual sequences that were accidently joined in the process somewhere.The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on on factors including experimental procedures and sample complexity.  DADA2 identifies likely chimeras by aligning each sequence with those that were recovered in greater abundance and then seeing if there are any lower-abundance sequences that can be made exactly by mixing left and right portions of two of the more-abundant ones. 
> ![image](https://github.com/user-attachments/assets/210b6caf-5af1-4a86-bbef-f58bfd62d880)
```R
# Removing chimeras -------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Identified 17 bimeras out of 2521 input sequences.

    # though we only lost 17 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.9931372 # in this case we barely lost any in terms of abundance

    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together
```
Let's take a look at what has happened to our input for each sample
```R
# Count Overview ----------------------------------------------------------

    # set a little function
getN <- function(x) sum(getUniques(x))

    # making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
               filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
               dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
               nonchim=rowSums(seqtab.nochim),
               final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab

# write this table to a file to save
write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)
 
```
## 🧪 Step 9: Assign taxonomy (finally the cool stuff!)
We are going to be using the alternative classification method with the DECIPHER package because it's performance and speed are reported to be better than DADA2's standard workflow.
However, this process takes ~20-30 min. If we are running low on time, we will skip this. 

```R
# Assign Taxonomy ---------------------------------------------------------

# downloading DECIPHER-formatted SILVA v138 reference
download.file(url="https://www2.decipher.codes/data/Downloads/TrainingSets/GTDB_r226-mod_April2025.RData", destfile="GTDB_r226-mod_April2025.RData")

# loading reference taxonomy object
load("GTDB_r226-mod_April2025.RData")

# creating DNAStringSet object of our ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

# and finally classifying. This will take about ~20 min. 
tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

```

If we did not have time, then download the classification/taxonomy table that I ran ahead of time for you using the Terminal in R. Go to Tools > Terminal > New Terminal. It will open next to your R console.

```bash
wget https://raw.githubusercontent.com/morgansobol/MicrobialGenomics-TXST-2025/main/data/03_16S_amplicon/tax_info.RData
```
If `wget` does not work, try curl instead:
```bash
curl -L -O https://raw.githubusercontent.com/morgansobol/MicrobialGenomics-TXST-2025/main/data/03_16S_amplicon/tax_info.RData
```

To load the object into R, do this:
```R
load("tax_info.RData") 
```

Ok, let's convert the output object of class "Taxa" to a taxonomic matrix analogous to the output if we had used the `assignTaxonomy` function instead. We will also generate a fasta file of our final ASV sequences and the count table. 

```R
# Extract the Output ------------------------------------------------------

  # giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
    asv_headers[i] <- paste(">ASV", i, sep="_")
}

    # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

    # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

    # tax table:
    # creating table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(tax_info, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
    }))
colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)
```
## 🧪 Bonus Step: Remove contaminants
Good science is all about good controls. Typically, sequencing projects should include a "blank" control, which is sample that never had DNA in it but went through the entire extraction and sequencing process. While we always hope the kits we use for such things are sterile, most of the times they actually are not 😬

We will use the package `decontam` in "prevalence mode to help us. In this method, the prevalence (presence/absence across samples) of each sequence feature in true positive samples is compared to the prevalence in negative controls to identify contaminants. here: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html 

So for `decontam` to work on our data, we need to provide it with our count table, currently stored in the “asv_tab” variable, and we need to give it a logical vector that tells it which samples are “blanks”. Here is making that vector and running the program:
```R
# Decontamination ---------------------------------------------------------

colnames(asv_tab) # our blanks are the first 4 of 20 samples in this case (B1, B2, B3, B4).
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 16))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 6 as contaminants

    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together

    # getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

asv_tax[row.names(asv_tax) %in% contam_asvs, ]
```
Since there were only six here, I wanted to peek at them. And not surprisingly, they are all things that are commonly contaminants, but of course not exclusively (e.g. Burkholderia, Pseudomonas). We can see this by looking at their taxonomic designations in our tax table. 

Let's also extract the contam sequences from the fasta files and blast against NCBI (a different database than GTDB) so we can possibly get an idea of what the NAs are. 
To do this, we can use the Terminal in R. 
```bash
grep -w -A1 "^>ASV_104\|^>ASV_219\|^>ASV_230\|^>ASV_274\|^>ASV_285\|^>ASV_623" ASVs.fa
```
Now paste into NCBI Blast https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome 

What did you find differently? 

And now, here is one way to remove the contaminants from our 3 primary outputs and create new files

```R
 # making new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

    # making new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

    # making new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

    ## and now writing them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
```

Now that concludes the end of the processing steps for 16S amplicon data in Dada2. Next exercise we will analyze our data and explore differences between samples statistically! 

## 📝 Assignment due before next lab on Canvas

Using the ASV_counts, ASV_taxonomy tables, and ASV.fa files, answer the following questions: 

1. What are the 3 most abundant taxa in the sample you chose?
2. Pick one of the three and tell me something about them (metabolism, lifestyle, other environments they are found, etc.). Provide references.
3. Blast the ASV you chose on NCBI. What were the top hits? Did they match?

