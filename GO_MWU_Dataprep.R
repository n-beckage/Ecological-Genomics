#### Go_MWU Data Prep ####
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
