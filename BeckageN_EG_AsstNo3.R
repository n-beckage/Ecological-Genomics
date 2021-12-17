library(ggplot2)
library(gridExtra)
library(GenomicRanges)
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)

setwd("C:/Users/nbeck/OneDrive/Desktop/Fall 2021/PBIO 395 Ecological Genomics/Homework/Assignment 3")

########### QUESTION 1 #############
# Get the list of admixed individuals:
Admixed <- read.table("Admixed.Inds",header=F)

# Get the meta data:
meta <- read.table("Combined_Transect_Sampling_Data_2020.txt", sep="\t",header=T)

# merge them together:
meta_admx <- merge(meta, Admixed, by.x="ID", by.y="V1")
str(meta_admx)  

# Read in the Admixture coefficients for KBals that we made from the K=5 file:
KBals <- read.table("Admixed_KBals", sep="\t", header=F)
names(KBals) = c("ID","KBals")

# Second merge:
meta_admx_KBals <- merge(meta_admx,KBals,by="ID")

# Bring in clim data:
clim <- read.table("climDat.txt",sep="\t",header=T)

# Merge clim data with meta and KBals:
meta_admx_KBals_clim <- merge(meta_admx_KBals,clim,by="ID")

#  Final freeze date
plotmean_finalFreeze <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=mean_finalFreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Average Date of Last Freeze") 

plotmean_finalFreeze

# Mean number of chilling degree days
plotmed_DD0 <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=med_DD0, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Average # of Annual Chilling Degree Days") 

plotmed_DD0

# Mean number of warming degree days
plotmean_cGDDfreeze <- ggplot(meta_admx_KBals_clim,aes(x=KBals,y=mean_cGDDfreeze, color=Transect.x)) +
  geom_point(size=2) + 
  xlab("Proportion P. balsamifera ancestry") +
  ylab("Average # of Annual Warming Degree Days") 

plotmean_cGDDfreeze

grid.arrange(plotmean_finalFreeze, plotmed_DD0, plotmean_cGDDfreeze, nrow = 3)

# linear models testing trait ~ genome-wide admixture association
summary(lm(mean_finalFreeze~KBals + Transect.x, data=meta_admx_KBals_clim))

summary(lm(med_DD0~KBals + Transect.x, data=meta_admx_KBals_clim))

summary(lm(mean_cGDDfreeze~KBals + Transect.x, data=meta_admx_KBals_clim))

######  Bring in Association results from Plink   ######
mean_finalFreeze <- read.table("plink2.mean_finalFreeze.glm.linear",skip=1,sep="\t",header=F)
names(mean_finalFreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
mean_finalFreeze2 <- mean_finalFreeze[which(mean_finalFreeze$TEST=="ADD"),]

# Define association outliers as the upper 0.1% of p-values

# Read in list of positions
snps <- read.table("Chr15.kept.sites",sep="\t", header=T)

#########  Average date of final freeze  #########
mean_finalFreeze2 <- cbind(snps, mean_finalFreeze2[,-c(1:2)])
mean_finalFreeze2$outlier = ifelse(mean_finalFreeze2$P<quantile(mean_finalFreeze2$P,0.01),2,1)
sum(mean_finalFreeze2$outlier==2)

p1 <- ggplot(mean_finalFreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=mean_finalFreeze2$outlier, color=mean_finalFreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Average date of last freeze")

p1

####### Average # of annual chilling degree days  #########
med_DD0 <- read.table("plink2.med_DD0.glm.linear",skip=1,sep="\t",header=F)
names(med_DD0) = c("CHROM",  "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
med_DD02 <- med_DD0[which(med_DD0$TEST=="ADD"),]
med_DD02 <- cbind(snps, med_DD02[,-c(1,2)])
med_DD02$outlier = ifelse(med_DD02$P<quantile(med_DD02$P,0.01),2,1)

p2 <- ggplot(med_DD02,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=med_DD02$outlier, color=med_DD02$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Average # of annual chilling degree days")

p2

#########  Average # of annual warming degree days  #########
mean_cGDDfreeze <- read.table("plink2.mean_cGDDfreeze.glm.linear",skip=1,sep="\t",header=F)
names(mean_cGDDfreeze) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
mean_cGDDfreeze <- mean_cGDDfreeze[which(mean_cGDDfreeze$TEST=="ADD"),]
mean_cGDDfreeze2 <- cbind(snps, mean_cGDDfreeze[,-c(1,2)])
mean_cGDDfreeze2$outlier = ifelse(mean_cGDDfreeze2$P<quantile(mean_cGDDfreeze2$P,0.01),2,1)

p3 <- ggplot(mean_cGDDfreeze2,aes(x=POS,y=-1*log10(P))) +
  geom_point(size=mean_cGDDfreeze2$outlier, color=mean_cGDDfreeze2$outlier) + 
  xlab("Position (bp) along chromosome") +
  ylab("-log10 P-value") +
  ggtitle("Average # of growing degree days")

p3

grid.arrange(p1, p2, p3, nrow = 3)

#########  Bud flush  #########
budflush <- read.table("plink2.FLUSH.glm.linear",skip=1,sep="\t",header=F)
names(budflush) = c("CHROM",    "POS",  "ID",   "REF",  "ALT",  "A1",   "TEST", "OBS_CT",   "BETA", "SE",   "T_STAT",   "P")
budflush <- budflush[which(budflush$TEST=="ADD"),]
budflush2 <- cbind(snps, budflush[,-c(1,2)])
budflush2$outlier = ifelse(budflush2$P<quantile(budflush2$P,0.01),2,1)

############ QUESTION 2 ###############
# Determining the ranges of 'peaks" for climate and bud flush variables
CHR="Chr15"
budflushGR <- GRanges(CHR,IRanges(budflush2$POS-2.5e4,budflush2$POS+2.5e4),POS=budflush2$POS, P=budflush2$P, outlier=budflush2$outlier)
mean_finalFreezeGR <- GRanges(CHR,IRanges(mean_finalFreeze2$POS-2.5e4,mean_finalFreeze2$POS+2.5e4),POS=mean_finalFreeze2$POS, P=mean_finalFreeze2$P, outlier=mean_finalFreeze2$outlier)
med_DD0GR <- GRanges(CHR,IRanges(med_DD02$POS-2.5e4,med_DD02$POS+2.5e4),POS=med_DD02$POS, P=med_DD02$P, outlier=med_DD02$outlier)
mean_cGDDfreezeGR <- GRanges(CHR,IRanges(mean_cGDDfreeze2$POS-2.5e4,mean_cGDDfreeze2$POS+2.5e4),POS=mean_cGDDfreeze2$POS, P=mean_cGDDfreeze2$P, outlier=mean_cGDDfreeze2$outlier)

# Subsetting to capture the outlier region - bud flush
budflushGRout <- unlist(reduce(split(budflushGR, ~outlier)))
budflushGRout$outlier <- names(budflushGRout)
budflushGRCand <- subset(budflushGRout, outlier==2)

# final freeze date
mean_finalFreezeGRout <- unlist(reduce(split(mean_finalFreezeGR, ~outlier)))
mean_finalFreezeGRout$outlier <- names(mean_finalFreezeGRout)
mean_finalFreezeGRCand <- subset(mean_finalFreezeGRout, outlier==2)

# mean annual chilling days
med_DD0GRout <- unlist(reduce(split(med_DD0GR, ~outlier)))
med_DD0GRout$outlier <- names(med_DD0GRout)
med_DD0GRCand <- subset(med_DD0GRout, outlier==2)

# Mean growing degree days
mean_cGDDfreezeGRout <- unlist(reduce(split(mean_cGDDfreezeGR, ~outlier)))
mean_cGDDfreezeGRout$outlier <- names(mean_cGDDfreezeGRout)
mean_cGDDfreezeGRCand <- subset(mean_cGDDfreezeGRout, outlier==2)

# Print the candidate regions
budflushGRCand
mean_finalFreezeGRCand
med_DD0GRCand
mean_cGDDfreezeGRCand

#### Overlap - each climatic variable and bud flush
overlap_BF_fFReeze <- subsetByOverlaps(budflushGRCand, mean_finalFreezeGRCand)
overlap_BF_DD0 <- subsetByOverlaps(budflushGRCand, med_DD0GRCand)
overlap_BF_cGDD <- subsetByOverlaps(budflushGRCand, mean_cGDDfreezeGRCand)

# Print the overlapping regions
overlap_BF_fFReeze
length(overlap_BF_fFReeze)
# [1] 0 -> no overlapping regions

overlap_BF_DD0
length(overlap_BF_DD0)
# [1] 1 -> one overlapping region

overlap_BF_cGDD
length(overlap_BF_cGDD)
# [1] 1 ->  one overlapping region

### Which genes are in these overlapping regions?
# Import the GFF annotation file and make a transcript database
txdb <- makeTxDbFromGFF("Ptrichocarpa_533_v4.1.gene.gff3.gz", format="gff3")
txdb

# Subset the database for just your chromosome of interest
seqlevels(txdb) <- CHR # subset for just your chromosome

# Reduce the transcript database to just the non-redundant gene names, instead of multiple entries for all the variant transcript types per gene
genes <- unlist(reduce(transcriptsBy(txdb, by="gene"))) 
genes$geneID <- names(genes)

# Finding candidate genes
candGenes_DD0 <- subsetByOverlaps(genes, overlap_BF_DD0)
candGenes_cGDD <- subsetByOverlaps(genes, overlap_BF_cGDD)

write.table(candGenes_DD0$geneID, paste0("candGenes_DD0_",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
write.table(candGenes_cGDD$geneID, paste0("candGenes_cGDD_",CHR,".txt"), quote=F, col.names=F, row.names=F, sep=",")
