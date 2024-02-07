##### This script takes your sequencing data and explore the quality of the data

##### Set up working directory #####
workdir <- '/scratch/Users/jewe1055/exp_files/exp83/output/'
setwd(workdir)
getwd()

##### Set up the environment #####
#install.packages("remotes")
#remotes::install_github("lpantano/DEGreport")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#install.packages("ggfortify")

library("DESeq2")
#library("DEGreport")
library(dplyr)    # alternatively, this also loads %>%,
library(ggplot2)
library(tidyr)    # gather()
library(tidyverse) 
library(ggfortify) #PCA and autoplot
library(ggrepel) #geom_text_repel
library(limma) #plotMA

##### Set up directories #####
outdir <- '/scratch/Users/jewe1055/exp_files/exp83/output/'
hg38gtf <- '/scratch/Shares/dowell/genomes/hg38/hg38_refseq.gtf'

##### Import counts table (featureCounts) and metadata table #####
#read the count table
coveragetable <- read.csv('/scratch/Users/jewe1055/exp_files/exp83/output/featureCounts/featureCounts_counts.csv', header=TRUE, sep="\t")
head(coveragetable)

#read the metadata
metadata <- read.csv("/scratch/Users/jewe1055/exp_files/exp83/scripts/metadata.txt", header=TRUE, sep=",")
head(metadata)

#check the count and metadata to make sure they are in the same order, have the same sample and drop any label that isn't in both tables
countdat <- coveragetable %>% select(as.vector(metadata$label))
head(countdat)

#set up the deseq object
dds <- DESeqDataSetFromMatrix(countData = countdat, colData = metadata, design = ~treatment+genotype)

#normalization factor
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

##### Quality control of samples via PCA
count_matrix <- as.matrix(countdat)
rownames(count_matrix) <- paste(coveragetable$GeneID) #rows as genes

#Filter counts < 5
keep <- rowSums(count_matrix) > 5
count_matrix <- count_matrix[keep,]
dim(count_matrix)
count_pca <- prcomp(t(count_matrix), scale=TRUE)

# summary of the prcomp object
summary(count_pca)

## get the name of the sample (cell) with the highest pc1 value
rownames(count_pca$x)[order(count_pca$x[ ,1], decreasing=TRUE)[1]]

## get the name of the top 10 measurements (genes) that contribute
## most to pc1.
loading_scores <- count_pca$rotation[,1]
gene_scores <- abs(loading_scores) ## get the magnitudes
gene_score_ranked <- sort(gene_scores, decreasing=TRUE)
top_10_genes <- names(gene_score_ranked[1:10])

top_10_genes ## show the names of the top 10 genes

## get the scree information
count_pca.var <- count_pca$sdev^2
scree <- count_pca.var/sum(count_pca.var)
plot((scree[1:10]*100), main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

##### Library sizes bar plot
librarySizes <- colSums(countdat)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")
abline(h=median(as.matrix(librarySizes)), col="red")

##### Count distribution
# Get log2 counts
logcounts <- log2(countdat + 1)

# Check distributions of samples using boxplots
boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2)
# Let's add a blue horizontal line that corresponds to the median
abline(h=median(as.matrix(logcounts)), col="blue")

# How does the count change with rlog()?
rlogcounts <- DESeq2::vst(count_matrix)
boxplot(rlogcounts, 
        xlab="", 
        ylab="rlog(Counts)",
        las=2)
abline(h=median(as.matrix(rlogcounts)), col="blue")

##### MA plot to see the composition bias
logcounts <- log2(countdat + 1)

par(mfrow=c(1,2))
plotMA(logcounts, array = 1)
abline(h=0, col="red")
plotMA(logcounts, array = 2)
abline(h=0, col="red")

# normalize the counts to solve for composition bias
normalizedCounts <- counts(dds, normalized=TRUE) 
logNormalizedCounts <- log2(normalizedCounts + 1)

par(mfrow=c(1,2))
plotMA(logNormalizedCounts, array = 1)
abline(h=0,col="red")
plotMA(logNormalizedCounts, array = 2)
abline(h=0,col="red")

#### Run PCA
pcDat <- prcomp(t(rlogcounts))
# plot PCA
autoplot(pcDat,
         data = metadata, 
         colour="name", 
         shape="treatment",
         size=5) +  
  geom_text_repel(aes(x=PC1, y=PC2, label=label), box.padding = 0.8)

# Remove the outlier to see how the PCA would plot
rlogcounts_edit <- subset(rlogcounts, select = -c(DR_BSA_2,Ursula_BSA_2, Ursula_IFNB_2))
metadata_edit <- metadata[-c(31, 43, 44), ] #31 DR_BSA2, 43 Ursula_BSA2, 44 Ursula_IFNB2
pcDat2 <- prcomp(t(rlogcounts_edit))

# plot PCA
autoplot(pcDat2,
         data = metadata_edit, 
         colour="name", 
         shape="treatment",
         size=5) +  
  geom_text_repel(aes(x=PC1, y=PC2, label=label), box.padding = 0.8)
