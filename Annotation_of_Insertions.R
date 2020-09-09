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
Insertions_GFF$gene<-gsub("chr","",Insertions_GFF$gene)
Insertions_GFF$Insertionsites<-gsub("chr","",Insertions_GFF$Insertionsites)
GFF3$gene<-paste(GFF3$V1,GFF3$V4,GFF3$V5)
Insertions_GFFinsidemerged<-merge(Insertions_GFF,GFF3,by="gene")
Insertions_GFFinside_geneIDs<-Insertions_GFFinsidemerged[,c(1,16,25)]

Mu_single_Ids<-merge(MuSingle,Insertions_GFFinside_geneIDs,by="Insertionsites")
Mu_single_Ids<-unique(Mu_single_Ids[,c(2,3,4,5,6,7,9)])
colnames(Mu_single_Ids)<-c("Chromosome","Start","End","Sample","StartReads","EndReads","Ids")
Mu_single_Ids<-Mu_single_Ids[grep("ID",Mu_single_Ids$Ids,invert = F),]
Mu_single_Ids<-Mu_single_Ids[grep("ID=gene:",Mu_single_Ids$Ids,invert = T),]
Mu_single_Ids<-Mu_single_Ids[grep("chromosome:",Mu_single_Ids$Ids,invert = T),]
Mu_single_Ids<-unique(Mu_single_Ids)
Mu_single_Ids$Ids<-gsub(";Parent=.*","",Mu_single_Ids$Ids)
Mu_single_Ids$Ids<-gsub("ID=","",Mu_single_Ids$Ids)

#create good gene id column
GeneID<-gsub("_T0.*","",Mu_single_Ids$Ids)
GeneID<-gsub("_T1.*","",GeneID)
GeneID<-gsub("_T2.*","",GeneID)
GeneID<-gsub("CDS:","",GeneID)
GeneID<-gsub("transcript:","",GeneID)
GeneID<-gsub("_P.*","",GeneID)
Mu_single_Ids$GeneID<-GeneID
Mu_single_Ids<-unique(Mu_single_Ids)

#create good transcript id column
TranscriptID <- gsub("CDS:","", Mu_single_Ids$Ids)
TranscriptID <- gsub("transcript:","", Mu_single_Ids$Ids)
Mu_single_Ids$TranscriptID<-TranscriptID
Mu_single_Ids<-unique(Mu_single_Ids)


#create two tables - one with gene IDs and one with TranscriptIds
#with unique we can get of the multiple repeating lines in the gene id dataframe
#here we remove all the lines with unnecessary content that are essentially doubles
#also remove all CDS lines from the Mu_single_TranscriptIds file - these are unnecessary doubles of the transcripts we don't need
Mu_single_GeneIds <- unique(Mu_single_Ids[, c("Chromosome", "Start", "End", "Sample", "StartReads", "EndReads", "GeneID")])
Mu_single_GeneIds <- Mu_single_GeneIds[grep("contig", Mu_single_GeneIds$GeneID, invert = TRUE), ]

Mu_single_TranscriptIds <- Mu_single_Ids[, c("Chromosome", "Start", "End", "Sample", "StartReads", "EndReads", "TranscriptID")]
Mu_single_TranscriptIds <- Mu_single_TranscriptIds[grep("CDS", Mu_single_TranscriptIds$TranscriptID, invert = TRUE), ]
Mu_single_TranscriptIds <- Mu_single_TranscriptIds[grep("contig", Mu_single_TranscriptIds$TranscriptID, invert = TRUE), ]

##
#add gene or transcript length to the dataframes
#gene length via gff3 file
#transcript length can be found in the gtf file
short_GFF3 <- GFF3[, c("V3", "V4", "V5", "V9")] 
#short_GFF3 <- short_GFF3[,-1]

#only keep rows with "gene" in column V3
short_GFF3_genes <- short_GFF3[grep("gene", short_GFF3$V3), ]

#delete everything but the gene name from column V9
short_GFF3_genes$V9 <- gsub("ID=gene:", "", short_GFF3_genes$V9)
short_GFF3_genes$V9 <- gsub("\\;.*", "", short_GFF3_genes$V9)

#compute gene length
Gene_length <- (short_GFF3_genes$V5 - short_GFF3_genes$V4 + 1)
short_GFF3_genes$Gene_length <- Gene_length

#write dataframe just with gene id and gene length and then merge it into the the Mu_single_gene_id dataframe
short_GFF3_genes <- short_GFF3_genes[, c("V9", "Gene_length")]
names(short_GFF3_genes) <- c("GeneID", "Gene_length")
Mu_single_GeneIds_gene_lengths <- merge(Mu_single_GeneIds, short_GFF3_genes, by="GeneID", all.x=TRUE)

###
#transcript length
#read in gtf file
MY_GTF <- list.files(getwd(), pattern="\\.gtf$") 
GTF <- read.delim(MY_GTF, header=FALSE, comment.char="#")


short_GTF <- GTF[, c("V3", "V4", "V5", "V9")] 

#only keep rows with "transcript" in column V3
short_GTF_transcripts <- short_GTF[grep("transcript", short_GTF$V3), ]

#delete everything but the transcript name from column V9
short_GTF_transcripts$V9 <- gsub("gene_id.*transcript_id ", "", short_GTF_transcripts$V9)
short_GTF_transcripts$V9 <- gsub("\\;.*", "", short_GTF_transcripts$V9)

#compute Transcript length
Transcript_length <- (short_GTF_transcripts$V5 - short_GTF_transcripts$V4 + 1)
short_GTF_transcripts$Transcript_length <- Transcript_length

#write dataframe just with gene id and gene length and then merge it into the the Mu_single_gene_id dataframe
short_GTF_transcripts <- short_GTF_transcripts[, c("V9", "Transcript_length")]
names(short_GTF_transcripts) <- c("TranscriptID", "Transcript_length")
Mu_single_TranscriptIds_transcript_lengths <- merge(Mu_single_TranscriptIds, short_GTF_transcripts, by="TranscriptID", all.x=TRUE)

####################
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

require(reshape2)
df <- melt(data.frame(x, y))
colnames(df) <- c("Row", "Column")

df_list <- as.list(df)

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
