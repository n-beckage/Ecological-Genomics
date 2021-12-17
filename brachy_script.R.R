library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn) 

setwd("C:/Users/nbeck/OneDrive/Desktop/Fall 2021/PBIO 395 Ecological Genomics/Group project")


#### Import the counts matrix and the metadata ####

# This is the counts matrix - it tells us how many reads get mapped to each transcript per sample.
countsTable <- read.table("brachy_countsMatrix_NO_HJ95.txt", header = T, row.names = 1) 
head(countsTable)
dim(countsTable)
# [1] 32415     8 -> 32,415 transcripts (genes) and 8 samples

# This step rounds all of the count values in our matrix so that we have whole numbers, since DESeq doesn't work with decimals. 
countsTable <- round(countsTable) 
head(countsTable)

# Then we import the metadata table, which tells us what treatment/conditions each sample received.
conds <- read.delim("brachy_samples.txt", header = T, stringsAsFactors = T, row.names = 1)
head(conds)

conds <- conds %>% 
  filter(Index!='HJ95s_24_S4')

# This histogram shows the distribution of the average number of reads each gene is mapped to. The large right hand tail indicates that there are very few genes that get expressed very often (average>1000); these are likely housekeeping genes.
hist(apply(countsTable,1,mean),xlim=c(0,5000),breaks=10000)

# Zooming in on the peak of this histogram, it appears that many genes are expressed 0-5 times. We will filter out genes that are expressed less than 10 times, as these are likely biologically insignificant
hist(apply(countsTable,1,mean),xlim=c(0,100),ylim=c(0,10000),breaks=100000)

# Median is a more informative data point than mean due to the large right hand tail
median(rowSums(countsTable))
# [1] 656

#### DESeq2 ####

# Create the DESeq data set
dds <- DESeqDataSetFromMatrix(countData = countsTable, colData = conds, design = ~Genotype)

dim(dds)
# We still all of our 32,415 transcripts, but now we have to filter

# Filter out transcripts with less than 10 reads per sample on average
dds <- dds[rowSums(counts(dds))>80]
dim(dds)
# brings us down to 21,897 transcripts

# Relevel such that WT genotype becomes the reference that the other genotypes are being compared too
#dds$Genotype <- relevel(dds$Genotype, ref = "Wild_Type")

# Run the model
dds <- DESeq(dds)

# Do the comparison of HDAC vs WT
HDAC_WT <- results(dds, contrast = c("Genotype","HDAC19","Wild_Type"),alpha= 0.05)
summary(HDAC_WT)
# Comparison of GF14H vs WT
GF14H_WT <- results(dds, contrast = c("Genotype","GF14H","Wild_Type"),alpha= 0.05)
summary(GF14H_WT)
# Might want to consider doing a third analysis for GF14H vs HDAC

# These lines will list, order, and summarize results (with p<0.05) from specific contrasts
GF14H_WT <- GF14H_WT[order(GF14H_WT$padj),]
HDAC_WT <- HDAC_WT[order(HDAC_WT$padj),]

# No N/A's!
# GF14H_WT <- GF14H_WT[!is.na(GF14H_WT$padj),]
# HDAC_WT <- HDAC_WT[!is.na(HDAC_WT$padj),]
dim(GF14H_WT)
summary(HDAC_WT)
summary(GF14H_WT)
# Removing NA's brings both down to 21,866 transcripts

# Shows the top-most significant DEGs for genotype HDAC knockout
head(HDAC_WT)
### FOR HDAC19 VS WT
# baseMean log2FoldChange     lfcSE      stat      pvalue
# <numeric>      <numeric> <numeric> <numeric>   <numeric>
#   Bradi4g41560.1  1936.761      -10.42774  0.683287 -15.26116 1.38769e-52
# Bradi1g21466.1   619.964        7.00278  0.507260  13.80511 2.37410e-43
# Bradi1g09177.1   506.161       -3.92787  0.358310 -10.96219 5.80745e-28
# Bradi5g03752.2   234.823        3.01212  0.288256  10.44947 1.47348e-25
# Bradi3g38140.1   119.719       -6.22486  0.627110  -9.92627 3.19997e-23
# Bradi2g54400.2   250.992       -4.61138  0.481293  -9.58123 9.58942e-22

head(GF14H_WT)
### FOR GF14H VS WT
# baseMean log2FoldChange     lfcSE      stat      pvalue
# <numeric>      <numeric> <numeric> <numeric>   <numeric>
#   Bradi1g21466.1  619.9638       -6.49368  0.508683 -12.76568 2.54919e-37
# Bradi3g35840.1  161.6551       -6.40262  0.655403  -9.76897 1.52995e-22
# Bradi4g43050.3   97.5821        3.62555  0.464020   7.81335 5.56861e-15
# Bradi3g09120.1 1328.7870       -2.85273  0.370630  -7.69698 1.39321e-14
# Bradi5g09920.2  471.3165        2.21111  0.291238   7.59209 3.14786e-14
# Bradi3g38140.1  119.7188        4.40289  0.583014   7.55194 4.28810e-14


summary(GF14H_WT)
# 185 genes (0.85%) are deferentially expressed in the GF14H knockout compared to the wild type; 0.37% (81) genes down-regulated, 0.48% (104) up-regulated

### Now for HDAC19
summary(HDAC_WT)
# 389 (1.8%) genes are differentially expressed in HDAC19 knockout genotype compared to the wild type; 132 (0.6%) are up-regulated, 257 (1.2%) are down-regulated

#### Go MWU Data Prep ####
# This script is meant to be appended to the project script that runs DESeq, as the following code only requires the two result objects from DESeq (GF14H_WT and HDAC_WT) as a priori inputs

library(dplyr)

# ------------- GF14h

# The names of genes of interest (DEGs), as character values
GF_geneID <- data.frame(row.names(GF14H_WT[,]))
GF_LFC <- data.frame(GF14H_WT[,2])
GF_geneID <- cbind(GF_geneID,GF_LFC)
names(GF_geneID) <- c("gene","logFoldChange")
write.csv(GF_geneID,file="GF14H_logFoldChange.csv",row.names=F)


# ----------- Gene Universe
# The gene Universe aka background set for GO analysis, which will contain the names of all 21897 transcripts
# Only need to make the geneUniverse.tab file once since it will be the same for both GF14h and HDAC19
transcriptName <- GF_geneID %>%
  select(1)
names(transcriptName)[1] <- "transcriptName"

# Read in the annotation file (with GO terms)
annotations <- read.delim("Bdistachyon_556_v3.2.annotation_info.txt")

# Selecting only the columns we need,transcript names and GO categories
GoMap <- annotations %>% 
  select(transcriptName,GO)

transcriptName <- merge(transcriptName,GoMap,all.x=T)


a=""
for (x in c(0:length(transcriptName[,2]))) {
  transcriptName[x,2] <- gsub(" ",";",transcriptName[x,2])
  if(identical(a,transcriptName[x,2])){
    transcriptName[x,2] <- gsub("","unknown",transcriptName[x,2])
  }
}

write.table(transcriptName, file="geneUniverse.tab", sep="\t", row.names = F,col.names = F)

# ------------- HDAC19
# The names of genes of interest (DEGs), as character values
HD_geneID <- data.frame(row.names(HDAC_WT[,]))
HD_LFC <- data.frame(HDAC_WT[,2])
HD_geneID <- cbind(HD_geneID,HD_LFC)
names(HD_geneID) <- c("gene","logFoldChange")
write.csv(HD_geneID,file="HDAC19_logFoldChange.csv",row.names=F)
