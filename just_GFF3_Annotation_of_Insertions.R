#####################################################################
# Annotation of previously identified Mutator insertion sites       #
# MuWU                                                              #
# Tyll Stoecker & Lena Altrogge                                     #
# tyll.stoecker@uni-bonn.de                                         #
#                                                                   #
# 20.03.2020                                                        #
#                                                                   #
#####################################################################

## This script adds various annoations to the previously created MuSeq SLI- .csv table
## It will create one file with gene level annotation and one with transcript level annotation
## It will further calculate the gene/transcript length and append this to the tables
## Lastly, we use the stock excel sheet of the MuSeq experiment to assess which F2 family is represented by the Row/Col intersection

library(ChIPpeakAnno)
library(IRanges)
library(GenomicRanges)
library(readxl)
library(reshape2)

options(warn=-1)

## Read table with identified Mu insertions
MuSingle <- read.csv("MuSeq_table_final/SLI-MuSeq_FGS.csv", header=T)
MuSingle$Insertionsites <- paste(MuSingle$Chr, MuSingle$Start, MuSingle$End)


## Read GFF File
setwd("FGS/")
MY_gff3 <- list.files(getwd(), pattern="\\.gff3$") 
GFF3 <- read.delim(MY_gff3, header=FALSE, comment.char="#")

GFF3 <- GFF3[GFF3$V3=="gene",]

## Calculate Ranges for Mu insertion sites and genes in the GFF
Ranges<-GRanges(seqnames=MuSingle$Chr,ranges=IRanges(start=MuSingle$Start,end=MuSingle$End))
GFF<-GRanges(seqnames = GFF3$V1,ranges=IRanges(start=GFF3$V4,end=GFF3$V5))

## Indentification of Insertions that are in genes in the GFF ##
annotatedPeak<-annotatePeakInBatch(Ranges, AnnotationData=GFF,output="inside")
Insertions_GFF <- as.data.frame(annotatedPeak)
Insertions_GFF <- na.omit(Insertions_GFF)

## add IDs to table with identified Mu insertion sites and write csv file 

Insertions_GFF$gene<-paste(Insertions_GFF$seqnames,Insertions_GFF$start_position,Insertions_GFF$end_position)
Insertions_GFF$Insertionsites<-paste(Insertions_GFF$seqnames,Insertions_GFF$start,Insertions_GFF$end)
#Insertions_GFF$gene<-gsub("chr","",Insertions_GFF$gene)
#Insertions_GFF$Insertionsites<-gsub("chr","",Insertions_GFF$Insertionsites)
GFF3$gene<-paste(GFF3$V1,GFF3$V4,GFF3$V5)
Insertions_GFFinsidemerged<-merge(Insertions_GFF,GFF3,by="gene")
Insertions_GFFinside_geneIDs<-Insertions_GFFinsidemerged[,c(1,16,25)]

Mu_single_Ids<-merge(MuSingle,Insertions_GFFinside_geneIDs,by="Insertionsites")
Mu_single_Ids<-unique(Mu_single_Ids[,c(2,3,4,5,6,7,9)])
colnames(Mu_single_Ids)<-c("Chromosome","Start","End","Sample","StartReads","EndReads","Ids")

Mu_single_Ids$Ids<-gsub("ID=gene:","",Mu_single_Ids$Ids)
Mu_single_Ids$Ids<-gsub(";.*","",Mu_single_Ids$Ids)


colnames(Mu_single_Ids)<-c("Chromosome","Start","End","Sample","StartReads","EndReads","GeneID")

#add gene or transcript length to the dataframes
#gene length via gff3 file
#transcript length can be found in the gtf file
short_GFF3 <- GFF3[, c("V3", "V4", "V5", "V9")] 
short_GFF3 <- short_GFF3[-1, ]

#delete everything but the gene name from column V9
short_GFF3$V9 <- gsub("ID=gene:","", short_GFF3$V9)
short_GFF3$V9 <- gsub(";.*","", short_GFF3$V9)


#compute gene length
Gene_length <- (short_GFF3$V5 - short_GFF3$V4 + 1)
short_GFF3$Gene_length <- Gene_length

#write dataframe just with gene id and gene length and then merge it into the the Mu_single_gene_id dataframe
short_GFF3_genes <- short_GFF3[, c("V9", "Gene_length")]
names(short_GFF3_genes) <- c("GeneID", "Gene_length")
Mu_single_GeneIds_gene_lengths <- merge(Mu_single_Ids, short_GFF3_genes, by="GeneID", all.x=TRUE)


#######
#transcript level

GFF3 <- read.delim(MY_gff3, header=FALSE, comment.char="#")

GFF3 <- GFF3[GFF3$V3=="mRNA",]

## Calculate Ranges for Mu insertion sites and genes in the GFF
Ranges<-GRanges(seqnames=MuSingle$Chr,ranges=IRanges(start=MuSingle$Start,end=MuSingle$End))
GFF<-GRanges(seqnames = GFF3$V1,ranges=IRanges(start=GFF3$V4,end=GFF3$V5))

## Indentification of Insertions that are in genes in the GFF ##
annotatedPeak<-annotatePeakInBatch(Ranges, AnnotationData=GFF,output="inside")
Insertions_GFF <- as.data.frame(annotatedPeak)
Insertions_GFF <- na.omit(Insertions_GFF)

## add IDs to table with identified Mu insertion sites and write csv file 

Insertions_GFF$gene<-paste(Insertions_GFF$seqnames,Insertions_GFF$start_position,Insertions_GFF$end_position)
Insertions_GFF$Insertionsites<-paste(Insertions_GFF$seqnames,Insertions_GFF$start,Insertions_GFF$end)
GFF3$gene<-paste(GFF3$V1,GFF3$V4,GFF3$V5)
Insertions_GFFinsidemerged<-merge(Insertions_GFF,GFF3,by="gene")
Insertions_GFFinside_geneIDs<-Insertions_GFFinsidemerged[,c(1,16,25)]

Mu_single_Ids<-merge(MuSingle,Insertions_GFFinside_geneIDs,by="Insertionsites")
Mu_single_Ids<-unique(Mu_single_Ids[,c(2,3,4,5,6,7,9)])
colnames(Mu_single_Ids)<-c("Chromosome","Start","End","Sample","StartReads","EndReads","Ids")


Mu_single_Ids$Ids<-gsub("ID=transcript:","",Mu_single_Ids$Ids)
Mu_single_Ids$Ids<-gsub(";.*","",Mu_single_Ids$Ids)


colnames(Mu_single_Ids)<-c("Chromosome","Start","End","Sample","StartReads","EndReads","TranscriptID")

#add gene or transcript length to the dataframes
#gene length via gff3 file
short_GFF3 <- GFF3[, c("V3", "V4", "V5", "V9")] 
short_GFF3 <- short_GFF3[-1, ]

#delete everything but the gene name from column V9
short_GFF3$V9 <- gsub("ID=transcript:","", short_GFF3$V9)
short_GFF3$V9 <- gsub(";.*","", short_GFF3$V9)


#compute transcript length
Transcript_length <- (short_GFF3$V5 - short_GFF3$V4 + 1)
short_GFF3$Transcript_length <- Transcript_length


#write dataframe just with gene id and gene length and then merge it into the the Mu_single_gene_id dataframe
short_GFF3_transcripts <- short_GFF3[, c("V9", "Transcript_length")]
names(short_GFF3_transcripts) <- c("TranscriptID", "Transcript_length")
Mu_single_TranscriptIds_transcript_lengths <- merge(Mu_single_Ids, short_GFF3_transcripts, by="TranscriptID", all.x=TRUE)


##########################################


#add stock information
setwd("../stock_matrix")
MY_STOCK <- list.files(getwd(), pattern="\\.xlsx$")
stock <- as.data.frame(read_excel(MY_STOCK))
names(stock)[1] <- "X__1"
rownames(stock) <- stock$X__1
stock <- stock[, -1]
stock <- as.matrix(stock)

#have to use the gene id otherwise run into problems with genes that are on top of each other (share the same loci)
#perform additional sorting as to circumvent wrong order
Mu_single_GeneIds_gene_lengths <- Mu_single_GeneIds_gene_lengths[order(Mu_single_GeneIds_gene_lengths[,1],
                                                                       Mu_single_GeneIds_gene_lengths[,2],
                                                                       Mu_single_GeneIds_gene_lengths[,3]),]

#create seperate dataframes - just with Cols and just with Rows
just_col <- Mu_single_GeneIds_gene_lengths[grep("Col_.*", Mu_single_GeneIds_gene_lengths$Sample), ]
just_row <- Mu_single_GeneIds_gene_lengths[grep("Row_.*", Mu_single_GeneIds_gene_lengths$Sample), ]

#merge the two based on everything but the sample and gene length column
test <- merge(just_row, just_col, by=c("GeneID", "Chromosome",  "Start",  "End"))

#add quotations
x = paste0('"', test$Sample.x, '"')
y = paste0('"', test$Sample.y, '"')

#install.packages("reshape2")
require(reshape2)
df <- melt(data.frame(x, y))
colnames(df) <- c("Row", "Column")

df_list <- as.list(df)

#this code has one drawback:
#if for some reason a row or column is never used it won't be part of the df_list
#this becomes obvious once the stock matrix and annotation is compared after MuWU is done
#> stock matrix annotation won't be corrected and everything is shifted
#this shouldn't be the case; however IF this should happen;
#check the factor level with str(df_list)
#and then add the missing factor level (at the right position)
#example: Row_01 never used; then we can insert easily in the beginning with the handy factor(levels()) funtion
#after that rerun this step of the workflow
#df_list$Row <- factor(df_list$Row, levels = c("\"Row_01\"", levels(df_list$Row)))


full <- c()
for (i in 1:nrow(test))  {
  full[i] <- stock[df_list$Row[i], df_list$Column[i]]
}

just_col$stock <- full
just_row$stock <- full

Mu_single_GeneIds_gene_lengths_and_stock <- (rbind(just_row, just_col))

#sort the dataframe so the two pairs are next to each other
Mu_single_GeneIds_gene_lengths_and_stock <- Mu_single_GeneIds_gene_lengths_and_stock[order(Mu_single_GeneIds_gene_lengths_and_stock[,1],
                                                                                           Mu_single_GeneIds_gene_lengths_and_stock[,2],
                                                                                           Mu_single_GeneIds_gene_lengths_and_stock[,3]),]


###########################################

#lastly add stock information to transcript level file
#for this we need a set of the complete geneId table but with only unique Sample, Start, End entries
#otherwise we get manymultiple matches and the table gets bigger
#this is still not perfect and we get some redundancies - but these are totally fine, since the gene level analysis is the focus
#and can always be used to check for each transcript
gene_merge_basis <- unique(Mu_single_GeneIds_gene_lengths_and_stock[, c("Start", "End", "Sample", "StartReads", "EndReads", "stock")])

Mu_single_TranscriptIds_transcript_lengths_and_stock <- merge(Mu_single_TranscriptIds_transcript_lengths, gene_merge_basis, 
                                                              by=c("Start", "End", "Sample", "StartReads", "EndReads"))


##
setwd("../MuSeq_table_final")
write.csv(Mu_single_GeneIds_gene_lengths_and_stock, "Mu_single_GeneIds_gene_lengths_and_stock.csv", row.names=F)
write.csv(Mu_single_TranscriptIds_transcript_lengths_and_stock, "Mu_single_TranscriptIds_transcript_lengths_and_stock.csv", row.names=F)
