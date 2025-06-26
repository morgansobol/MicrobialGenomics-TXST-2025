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
