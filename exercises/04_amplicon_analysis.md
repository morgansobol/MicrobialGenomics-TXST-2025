# Week XX: 16S Amplicon analysis - Part 1

Welcome to Week XX! 
This week, you will finally start analyzing data, specially performing 16S amplicon analysis. From here on out, unless noted, we are working in R, not at the Unix-like command line. If you find you need some refresher in R basics check out this page: https://astrobiomike.github.io/R/.

We will be following Dada2's tutortial + Mike Lee's tutorial on Dada2, using Mike's data (with some modifications) below. Thanks Mike!

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

Write a table/data matrix as a tab-delimited file
```R
write.table(object, "filename.tsv", sep="\t")
```
## 🧪 Step 1: Setting up the working environment
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
## 🧪 Step 2: More quality trimming/filtering
We did a filtering step above with cutadapt (where we eliminated reads that had imperfect or missing primers and those that were shorter than 215 bps or longer than 285), but in DADA2 we’ll implement a trimming step as well (where we trim reads down based on some quality threshold rather than throwing the read away). 

Since we’re potentially shortening reads further, we’re again going to include another minimum-length filtering component. We can also take advantage of a handy quality plotting function that DADA2 provides to visualize how you’re reads are doing, `plotQualityProfile()`. By running that on our variables that hold all of our forward and reverse read filenames, we can easily generate plots for all samples or for a subset of them. So let’s take a peak at that to help decide our trimming lengths:

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

>[!TIP]
> Zymo developed a tool called Figaro (I have not used it myself yet) that can help you choose DADA2 trimming  parameters: https://github.com/Zymo-Research/figaro#figaro

In DADA2, this quality-filtering step is done with the `filterAndTrim()` function:
```R
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                reverse_reads, filtered_reverse_reads, maxEE=c(2,2),
                rm.phix=TRUE, minLen=175, truncLen=c(250,200))
```
Here, the first and third arguments (“forward_reads” and “reverse_reads”) are the variables holding our input files, which are our primer-trimmed output fastq files from cutadapt. The second and fourth are the variables holding the file names of the output forward and reverse seqs from this function. And then we have a few parameters explicitly specified:
- `maxEE` is the quality filtering threshold being applied based on the expected errors and in this case we are saying we want to throw the read away if it is likely to have more than 2 erroneous base calls (we are specifying for both the forward and reverse reads separately).
- `rm.phix` removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring. And minLen is setting the minimum length reads we want to keep after trimming.
- As mentioned above, the trimming occurring beyond what we set with truncLen is coming from a default setting, `truncQ`, which is set to 2 unless we specify otherwise, meaning it trims all bases after the first quality score of 2 it comes across in a read.
- There is also an additional filtering default parameter that is removing any sequences containing any Ns, `maxN`= 0 by default.
- Then we have our `truncLen` parameter setting the minimum size to trim the forward and reverse reads to in order to keep the quality scores roughly above 30 overall.

As mentioned, the output read files were named in those variables we made above (“filtered_forward_reads” and “filtered_reverse_reads”), so those files were created when we ran the function – we can see them if we run `list.files()` in R, or by checking in our working directory in the terminal:

Let's check the object we generated called *filtered_out*, containing how many reads went in and how many reads made it out. How do we view this object?

[insert pic]

Now plot the reads like before using the `plotQualityProfile` function, *but* make sure you do it for the filtered reads this time. How do we do that? 

[insert pic]

## 🧪 Step 3: Generating an error model of our data
The dada2 tool that we will use for this analysis uses a statistic approach that aims to predict if sequences are actually true biological sequences or are partly generated by errors that appear during sequencing.

Each sequencing run, even when all goes well, will have its own subtle variations to its error profile. This step tries to assess that for both the forward and reverse reads. It is one of the more computationally intensive steps of the workflow. We will use the `multithread=TRUE` to speed this up. 
```R
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)
```
Now we plot these results to see what errors we have
```R
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
```
The red line is what is expected based on the quality score, the black line represents the estimate, and the black dots represent the observed. Generally speaking, you want the observed (black dots) to track well with the estimated (black line). So things look good and we can move on!

## 🧪 Step 4: Inferring ASVs
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

## 🧪 Step 5: Merging forward and reverse reads
Now DADA2 merges the forward and reverse ASVs to reconstruct our full target amplicon requiring the overlapping region to be identical between the two. By default it requires that at least 12 bps overlap, but in our case the overlap should be much greater. If you remember above we trimmed the forward reads to 250 and the reverse to 200, and our primers were 515f–806r. After cutting off the primers we’re expecting a typical amplicon size of around 260 bases, so our typical overlap should be up around 190. However, that’s estimated based on E. coli 16S rRNA gene positions and very back-of-the-envelope-esque of course, so to allow for true biological variation and such I’m going ot set the minimum overlap for this dataset for 170. I’m also setting the trimOverhang option to TRUE in case any of our reads go passed their opposite primers (which I wouldn’t expect based on our trimming, but is possible due to the region and sequencing method).

```R
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
## 🧪 Step 6: Generating a count table
Now we can generate a count table with the makeSequenceTable() function. This is one of the main outputs from processing an amplicon dataset. You may have also heard this referred to as a biome table, or an OTU matrix.
```R
seqtab <- makeSequenceTable(merged_amplicons)
class(seqtab) # matrix
dim(seqtab) # 20 2521

    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together
```
The count table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. 

## 🧪 Step 7: Remove chimeras
Chimeras are separate, individual sequences that were accidently joined in the process somewhere.The frequency of chimeric sequences varies substantially from dataset to dataset, and depends on on factors including experimental procedures and sample complexity.  DADA2 identifies likely chimeras by aligning each sequence with those that were recovered in greater abundance and then seeing if there are any lower-abundance sequences that can be made exactly by mixing left and right portions of two of the more-abundant ones. 
> ![image](https://github.com/user-attachments/assets/210b6caf-5af1-4a86-bbef-f58bfd62d880)
```R
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T) # Identified 17 bimeras out of 2521 input sequences.

    # though we only lost 17 sequences, we don't know if they held a lot in terms of abundance, this is one quick way to look at that
sum(seqtab.nochim)/sum(seqtab) # 0.9931372 # in this case we barely lost any in terms of abundance

    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together
```
Let's take a look at what has happened to our input for each sample
```R
    # set a little function
getN <- function(x) sum(getUniques(x))

    # making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
               filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
               dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
               nonchim=rowSums(seqtab.nochim),
               final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

summary_tab
#       dada2_input filtered dada_f dada_r merged nonchim final_perc_reads_retained
# B1           1613     1498   1458   1466   1457    1457                      90.3
# B2            591      529    523    524    523     523                      88.5
# B3            503      457    450    451    450     450                      89.5
# B4            507      475    440    447    439     439                      86.6
# BW1          2294     2109   2066   2082   2054    2054                      89.5
# BW2          6017     5527   5134   5229   4716    4716                      78.4
# R10         11258    10354   9658   9819   9009    8847                      78.6
# R11BF        8627     8028   7544   7640   7150    6960                      80.7
# R11          8927     8138   7279   7511   6694    6577                      73.7
# R12         15681    14423  12420  12932  10714   10649                      67.9
# R1A         12108    10906   9584   9897   8559    8535                      70.5
# R1B         16091    14672  12937  13389  11202   11158                      69.3
# R2          17196    15660  14039  14498  12494   12436                      72.3
# R3          17494    15950  14210  14662  12503   12444                      71.1
# R4          18967    17324  16241  16501  14816   14750                      77.8
# R5          18209    16728  14800  15332  12905   12818                      70.4
# R6          14600    13338  11934  12311  10459   10448                      71.6
# R7           8003     7331   6515   6726   5630    5618                      70.2
# R8          12211    11192  10286  10513   9530    9454                      77.4
# R9           8600     7853   7215   7390   6740    6695                      77.8

    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together

write.table(summary_tab, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)
   ## write this table to a file to save
```
## 🧪 Step 8: Assign taxonomy (finally the cool stuff!)
We are going to be using the alternative classification method with the DECIPHER package because it's performance and speed are reported to be better than DADA2's standard workflow.
```R
## load package
library(DECIPHER); packageVersion("DECIPHER")

## downloading DECIPHER-formatted SILVA v138 reference
download.file(url="https://www2.decipher.codes/data/Downloads/TrainingSets/SILVA_SSU_r138_2_2024.RData", destfile="SILVA_SSU_r138_2_2024.RData")

## loading reference taxonomy object
load("SILVA_SSU_r138_2_2024.RData")

## creating DNAStringSet object of our ASVs
# dna <- DNAStringSet(getSequences(seqtab.nochim))

## and classifying
# tax_info <- IdTaxa(test=dna, trainingSet=trainingSet, strand="both", processors=NULL)

```

Ok, let's convert the output object of class "Taxa" to a taxonomic matrix analogous to the output if we had used the `assignTaxonomy` function instead. We will also generate a fasta file of our final ASV sequences and the count table. 

```R
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

We will use the package `decontam` to help us. 
```R
library(decontam)
packageVersion("decontam")
#
```
For `decontam` to work on our data, we need to provide it with our count table, currently stored in the “asv_tab” variable, and we need to give it a logical vector that tells it which samples are “blanks”. Here is making that vector and running the program:
```R
colnames(asv_tab) # our blanks are the first 4 of 20 samples in this case
vector_for_decontam <- c(rep(TRUE, 4), rep(FALSE, 16))

contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

table(contam_df$contaminant) # identified 6 as contaminants

    ## don't worry if the numbers vary a little, this might happen due to different versions being used 
    ## from when this was initially put together

    # getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])
```
Since there were only six here, I wanted to peek at them. And not surprisingly, they are all things that are commonly contaminants, but of course not exclusively (e.g. Burkholderia, Escherichia, Pseudomonas, Corynebacterium). We can see this by looking at their taxonomic designations in our tax table:
```R
asv_tax[row.names(asv_tax) %in% contam_asvs, ]

#         Kingdom    Phylum           Class                 Order                   Family               Genus                                       
# ASV_104 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Betaproteobacteriales" "Burkholderiaceae"   "Burkholderia-Caballeronia-Paraburkholderia"
# ASV_219 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Enterobacteriales"     "Enterobacteriaceae" "Escherichia/Shigella"                      
# ASV_230 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pseudomonadales"       "Pseudomonadaceae"   "Pseudomonas"                               
# ASV_274 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Pseudomonadales"       "Pseudomonadaceae"   "Pseudomonas"                               
# ASV_285 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Betaproteobacteriales" "Burkholderiaceae"   "Tepidimonas"                               
# ASV_623 "Bacteria" "Actinobacteria" "Actinobacteria"      "Corynebacteriales"     "Corynebacteriaceae" "Corynebacterium_1"

    ## again don't worry if the numbers vary a little, this might happen due to different versions being used from when this was initially put together
```
And now, here is one way to remove them from our 3 primary outputs and create new files

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

## 📝 Assignment due next class on Canvas
Perform a top-level exploration of one of the sample outputs and answer these three questions. 
1. What are the dominant taxa in the sample you chose?
2. Tell me something about what this taxa is known for.
3. What % of your sample is unclassified? (I.e. not assignment to the genus level). 
