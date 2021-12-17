#### Setup ####
# First step is to load all of the necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(BiocManager)
library(DESeq2)
library(vsn)
library(pheatmap)
library(eulerr)

# Then set the working directory to the location of the data files.
setwd("C:/Users/nbeck/Desktop/College & Career/Classes/4thYearSchoolWork/Fall 2021/PBIO 395 Ecological Genomics/RNAseq/tonsa_play/Assignment 2")

#### Import the counts matrix and the metadata ####

# This is the counts matrix - it tells us how many reads get mapped to each transcript per sample.
countsTable_F3 <- read.table("DE_counts_F3.txt", header = T, row.names = 1) 
head(countsTable_F3)
dim(countsTable_F3)
# [1] 25279    16 -> 25,279 transcripts (genes) and 16 samples

# This step rounds all of the count values in our matrix so that we have whole numbers, since DESeq doesn't work with decimals. 
countsTable_F3 <- round(countsTable_F3) 
head(countsTable_F3)

# Then we import the metadata table, which tells us what treatment/conditions each sample received.
conds_F3 <- read.delim("RT_tonsa_F3_samples.txt", header = T, stringsAsFactors = T, row.names = 1)
head(conds_F3)

# This histogram shows the distribution of the average number of reads each gene is mapped to. The large right hand tail indicates that there are very few genes that get expressed very often (average>1000); these are likely housekeeping genes.
hist(apply(countsTable_F3,1,mean),xlim=c(0,5000),breaks=10000)

# Zooming in on the peak of this histogram, it appears that most genes are expressed 20-40 times on average; few are expressed less than 20 times on average. So, I will use 20 as my cutoff when I create the DESeq object below.
hist(apply(countsTable_F3,1,mean),xlim=c(0,100),ylim=c(0,3000),breaks=14000)

# Median is a more informative data point than mean due to the large right hand tail
median(rowSums(countsTable_F3))
# [1] 2144

#### Interaction effects ####

# This object quantifies the effect of copepod evolutionary line, the environment they were reciprocated in, and the interaction between these two factors.
dds <- DESeqDataSetFromMatrix(countData = countsTable_F3, colData = conds_F3,
                              design = ~ line + environment + line:environment)

# Now I filter out genes with too few reads. I decide to keep reads with an average > 20 reads per sample, as explained with the histogram above.
dds <- dds[rowSums(counts(dds))>320]

# As we can see we filtered out 744 genes, which is the quantity matching the leftmost bin on our histogram.
dim(dds)
# [1] 24536    16

# This line runs the DESeq model to test for differential gene expression due to environment + line effects (interaction) with our data set using the likelihood ratio test
dds <- DESeq(dds, test="LRT", reduced=~environment + line)

# These lines will list, order, and summarize results (with p<0.05) from specific contrasts
resInt <- results(dds, alpha = 0.05)
resInt <- resInt[order(resInt$padj),]

# Shows the top-most significant DEGs for interaction effects.
head(resInt)
# log2 fold change (MLE): linecombined.environmentHH 
# LRT p-value: '~ line + environment + line:environment' vs '~ environment + line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN142181_c0_g4    384.575        3.12233  0.396876   59.9355 9.80184e-15 2.12465e-10
# TRINITY_DN143012_c0_g4    468.423       -3.31494  0.476170   46.2058 1.06461e-11 1.15383e-07
# TRINITY_DN131723_c0_g1   1927.579        2.61189  0.386190   44.3215 2.78642e-11 2.01328e-07
# TRINITY_DN142181_c0_g18   365.047        3.12301  0.466624   43.1913 4.96401e-11 2.68999e-07
# TRINITY_DN145818_c5_g1    298.152        1.90159  0.297244   40.3943 2.07544e-10 8.99743e-07
# TRINITY_DN135177_c0_g1   5857.774        2.38746  0.395992   35.3690 2.72793e-09 9.85508e-06


# 0.93% (229) genes up-regulated, 0.21% (51) down-regulated as a result of interaction. 1.14% (280) total DEGs.
summary(resInt)
# out of 24536 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 229, 0.93%
# LFC < 0 (down)     : 51, 0.21%
# outliers [1]       : 5, 0.02%
# low counts [2]     : 2855, 12%
# (mean count < 31)

# No N/A's!
resInt <- resInt[!is.na(resInt$padj),]

# Isolating DEGs due to interaction
degsInt <- row.names(resInt[resInt$padj < 0.05,])

# Should match our sum from the summary above, and indeed it does.
length(degsInt)
# [1] 280

head(degsInt)

#### Environment effects ####

# This line sets up the design to explicitly look for differences in gene expression across environmental conditions.
dds <- DESeqDataSetFromMatrix(countData = countsTable_F3, colData = conds_F3, 
                              design = ~ line + environment)

# Same filter step as above
dds <- dds[rowSums(counts(dds))>320]

# Running the model with a likelihood ratio test ("LRT")
dds <- DESeq(dds, test="LRT", reduced=~line)

# These lines will list, order, and summarize results (with p<0.05) from specific contrasts
resEnv <- results(dds, alpha = 0.05)
resEnv <- resEnv[order(resEnv$padj),]

# Shows the top-most significant DEGs for environmental effect
head(resEnv)
# log2 fold change (MLE): environment HH vs AA 
# LRT p-value: '~ line + environment' vs '~ line' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN121599_c1_g1     728.891       -6.01655  0.471623  109.3964 1.32872e-25 3.19478e-21
# TRINITY_DN150588_c1_g2    1342.265        2.10049  0.271438   53.6328 2.41695e-13 2.90566e-09
# TRINITY_DN136932_c13_g13   133.132        2.10156  0.292811   47.8211 4.66949e-12 3.74244e-08
# TRINITY_DN146851_c0_g5     465.895        1.33883  0.193292   46.5271 9.03583e-12 5.43143e-08
# TRINITY_DN145745_c0_g5     130.619        1.46458  0.218487   43.6640 3.89869e-11 1.87480e-07
# TRINITY_DN83766_c0_g1      109.786       -1.56734  0.233268   43.2931 4.71240e-11 1.88842e-07

# 2.1% of genes (512) were up-regulated, 1.3% (318) down-regulated in response to hi-temp, hi-CO2 environment compared to ambient. 3.4% (830) differentially expressed in total.
summary(resEnv)
# out of 24536 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 512, 2.1%
# LFC < 0 (down)     : 318, 1.3%
# outliers [1]       : 16, 0.065%
# low counts [2]     : 476, 1.9%
# (mean count < 22)

# filters matrix to exclude N/A's
resEnv <- resEnv[!is.na(resEnv$padj),]

# Take out the genes that were differentially expressed across environmental treatments
degsEnv <- row.names(resEnv[resEnv$padj < 0.05,]) 

# The length of degsEnv should match the number of differentially expressed genes as indicated in the model summary above. Indeed it does.
length(degsEnv)
# [1] 830

#### Line effects ####

# Looking at effects of copepod line now; "line" goes second after tilda
dds <- DESeqDataSetFromMatrix(countData = countsTable_F3, colData = conds_F3, 
                              design = ~ environment + line)

# Can't forget to filter!
dds <- dds[rowSums(counts(dds))>320]

# Now run the model, same as before except we're comparing against environment (reduced=~environment)
dds <- DESeq(dds, test="LRT", reduced=~environment)

# List, order and summarize again
resLine <- results(dds, alpha = 0.05)
resLine <- resLine[order(resLine$padj),]

# Shows the top-most significant DEGs for line effect
head(resLine)
# log2 fold change (MLE): line combined vs ambient 
# LRT p-value: '~ environment + line' vs '~ environment' 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   TRINITY_DN132194_c0_g1  1230.730        1.53627  0.154432   94.8431 2.06087e-22 4.95516e-18
# TRINITY_DN121089_c0_g2   633.118        1.49854  0.163992   80.3309 3.16672e-19 3.80703e-15
# TRINITY_DN134798_c0_g1   153.246       -2.05474  0.224450   77.8584 1.10701e-18 8.87228e-15
# TRINITY_DN129890_c0_g4   154.223       -3.22656  0.341745   76.6075 2.08547e-18 1.25358e-14
# TRINITY_DN147342_c0_g4   151.670        1.86447  0.238133   58.4042 2.13440e-14 1.02639e-10
# TRINITY_DN134960_c1_g9  2872.433        1.93785  0.245326   57.2969 3.74746e-14 1.50173e-10

# 3.3% (805) genes up-regulated, 3.4% (832) down-regulated across evolutionary lines controlling for differences in environment and interaction. 6.7% (1637) genes DE total.
summary(resLine)
# out of 24536 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 805, 3.3%
# LFC < 0 (down)     : 832, 3.4%
# outliers [1]       : 16, 0.065%
# low counts [2]     : 476, 1.9%
# (mean count < 22)

# Filtering out N/A's
resLine <- resLine[!is.na(resLine$padj),]

# Take out DEGs significant in response to evolutionary lines
degsline <- row.names(resLine[resLine$padj < 0.05,])

# The length of degsLine should match the number of DEGs as indicated in the model summary above. Checks out.
length(degsline)
# [1] 1637

head(degsline)


#### VENN DIAGRAM ####

# Total numbers of differentially expressed genes by effect
length(degsEnv)  # 830
length(degsline)  # 1637
length(degsInt)  # 280

# Intersections ->  overlap (in Venn diagram speak) between groups
intEL <- length(intersect(degsEnv,degsline))  # 139
intEL
intEI <- length(intersect(degsEnv,degsInt))  # 14
intEI
intLI <- length(intersect(degsInt,degsline))  # 31
intLI

# Overlap of all three groups for venn diagram
intersectEL <- intersect(degsEnv,degsline)
intELI <- length(intersect(intersectEL,degsInt))
intELI
# [1] 7

# Number unique for each group -> these numbers will be in only in the unique circle of each respective group in the venn diagram
E <- 830-139-14-7
L <- 1637-139-31-7
I <- 280-14-31-7

# Creating the venn object
fit1 <- euler(c("Env" = E, "Line" = L, "Interaction" = I, "Env&Line" = intEL, "Env&Interaction" = intEI, "Line&Interaction" = intLI, "Env&Line&Interaction" = intELI))

# Always nice to finish a long R script with a pretty colorful plot courtesy of Wes Anderson
plot(fit1, quantities = TRUE, fills=wes_palette("GrandBudapest2"))
