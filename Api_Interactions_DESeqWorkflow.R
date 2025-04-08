#R script for Differential Expression Analysis workflow
#https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-count-matrix

#For Differential Expression

#-----------Load libraries

library(DESeq2)
library(tidyverse)
library(dplyr)
library(PCAtools)
library(org.Mm.eg.db)
library(data.table)
library(matrixStats)
library("RColorBrewer")
library(pheatmap)
library(openxlsx)

#From the DESEQ vignette - if we want to test weather a variable improves the model we can use the likelihood ratio test
#So for the scope of testing wether the term Sex:Group is relevant to the model and which genes are mostly affected by this we will use the LRT



#-----------Define Variables
#Comprehensive directory
Mum_dir='W:/01_GenEpi_Projects/Schulte_2025/4.Interactions_CTR_asREF'

wdir=Mum_dir
##Control group to use as reference level for the DESeq Object
CTRL_level='CTR'

## FDR cutoff for results() and results summary 
#padj threshold for statistical significance (for alpha=)
alpha_var =0.05
#Log2FoldChange threshold for the Wald test (for lfcThreshold=). Default altHypothesis tests for
#LFC > or < of absolute value specified by lfcThreshold=
lfcThreshold_var=0

#cutoffs for significance
#Fold change threshold in absolute value
FC_thr=0
#pvalue thresholds
pval_thr=0.05


#-----------Load data+ check data format

#Set working directory
setwd(paste(wdir))

#Input data: count matrix generated with Htseq-count where the column names correspond to sample names, in the same order
#as the row names in the cov matrix
count_matrix <- read.table('CountMatrix.txt', sep='\t', header=TRUE, row.names=1)

#Matrix containing the metadata for the samples, row names correspond to column names in count_matrix
#Contains the 'Group' variable to be used in the experimental design and that includes the reference level for comparisons
cov <- read.table('Groups.txt', sep='\t', header=TRUE, row.names=1)

#Subset dataset to use only samples in cov and in the same order
#Create subset of transformed dataset for the samples of interest 
library(dplyr)

#Create vector for subsetting the dataset for the samples of interest
Subset_v <- row.names(cov)

#Subset dataset of transformed data for heatmap
count_matrix<- subset(count_matrix, select = Subset_v)
head(count_matrix)

#Making sure that row names in cov data match the columns in the count matrix
all(rownames(cov) %in% colnames(count_matrix))
#are they in the same order?
all(rownames(cov) == colnames(count_matrix))

#Create GENEID data frame file for gene annotation
#Get the ENSEMBL IDs as vector
#library(org.Mm.eg.db)
ENSEMBLIDS <- count_matrix
ENSEMBLIDS$ENSEMBL <- rownames(count_matrix)
ENSEMBLIDS<-as.vector(ENSEMBLIDS$ENSEMBL)
#Use org.Mm.eg.db package to map Gene Symbol (or other) to the ENSEMBL IDs in the GENEID dataframe
GENEID <- as.data.frame (mapIds(org.Mm.eg.db, ENSEMBLIDS, keytype="ENSEMBL", column="SYMBOL", multiVals = "first"))
GENEID$ENSEMBL <- rownames(GENEID)
#Rename column to GeneSymbol (GENEID columns: 'GeneSymbol', 'ENSEMBL')
#library(data.table)
setnames (GENEID, paste(colnames(GENEID[1])), 'GeneSymbol')
colnames (GENEID[1]) <- 'GeneSymbol'



#-----------Construct the DESeq object 'dds' with experimental design
############ Define Model Design
############ The variable of interest goes at the END of the formula (right-hand side)
#this way the results function will by default pull the results in relation to that variable
#It is also possible to retrieve the log2 fold changes, p values and adjusted p values of variables 
#other than the last one in the design.With 'contrast' 

dds<-DESeqDataSetFromMatrix(countData = count_matrix,
                            colData=cov,
                            design= ~ Sex + Group + Sex:Group)
dds

#When we add an interaction variable then the 'Group' variable results are relative only for the reference level 
#so we will use this design only to perform a likelihood ratio test 

#NK: Next, we estimate the size factors (which is a representation of the differences in coverage 
#between replicates):

ndata <- estimateSizeFactors(dds)
sizeFactors(ndata)

#pre-filtering:remove rows with low gene counts
#According to Nils and Andreas it is preferable to keep the dataset intact
#for as long as possible, especially if we are not limited by power of the computer
#However, at the end we have many rows with empty counts and very low counts
#so I'd rather filter it as following the Biocondoctur vignette, also since as stated here
#it will improve visualisation

smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

dds

#------------Define control group
#set the factor level
#If we don't explicitly chose a reference level it will just select alphabetically and use that as reference
dds$Group <-relevel (dds$Group, ref = CTRL_level)
dds$Sex <- relevel (dds$Sex, ref = 'M')

#Changing the reference level does not give us additional genes, it only affects the log2foldchange (because we are changing the orders of numerator and denominator)
#dds$Group <-relevel (dds$Group, ref = CTRL_level)
#dds$Group <-relevel (dds$Group, ref='TGxKO')
#dds$Sex <- relevel (dds$Sex, ref = 'F')


#------------Run DESeq function to apply the LRT
#In this step the data is normalised by the size factor
#Express the reduced design where the variable to test has been removed 
dds <-DESeq(dds, test='LRT', reduced = ~Sex + Group )
resultsNames(dds)

  #"SexF.GroupCremTG"
  #"SexF.GroupHDAC2_KO"
  #"SexF.GroupTGxKO"

{
  res<-results(dds, name =  "SexF.GroupTGxKO")
  res
  summary(res)
  resultsNames(dds)

  ####Save annotated table
  #Create result data frame
  res_df<-as.data.frame(res)

  #Annotate with GeneSymbol column and rearrange columns
  res_df <- merge(res_df, GENEID, by=0)
  #replace missing GeneSymbol (NAs) with original ENSEMBL ID.
  res_df<-res_df %>%
   mutate(GeneSymbol = coalesce(GeneSymbol, ENSEMBL)) # works


  #remove additional columns and reorder
  lastcol<-(ncol(res_df))
  res_df<-res_df[,2:lastcol]  
  res_df <- res_df %>%
   dplyr::select(ENSEMBL, GeneSymbol, everything()) 

  res_df <-  res_df %>%
    arrange(padj)  # arrange in ascending order, NAs are put at the bottom of the file
  head(res_df)

####FOR HEATMAP: Create a subset dataset with only genes with padj <pval_thr sorted by decreasing |Shrunk_LFC|
#Keep only genes with padj <pval_thr
  res_df_sig<- res_df[res_df$padj <= pval_thr,]
#Remove rows where padj = missing 
  res_df_sig <- res_df_sig %>% drop_na(padj)
  res_df_sig <-  res_df_sig %>%
   arrange(padj)  # arrange in ascending order, NAs are put at the bottom of the file
  head(res_df_sig)
}


write.table(res_df,"DESeq2_Results_Interaction_Sex-Group_TGKOvsCTR_FvsM.txt",sep="\t", quote=F,row.names = FALSE)
write.xlsx(res_df,"DESeq2_Results_Interaction_Sex-Group_TGKOvsCTR_FvsM.xlsx")
write.xlsx(res_df_sig,"DESeq2_Results_padj_Interaction_Sex-Group_TGKOvsCTR_FvsM.xlsx")

write.table(res_df,"DESeq2_Results_Interaction_Sex-Group_CremTGvsCTR_FvsM.txt",sep="\t", quote=F,row.names = FALSE)
write.xlsx(res_df,"DESeq2_Results_Interaction_Sex-Group_CremTGvsCTR_FvsM.xlsx")
write.xlsx(res_df_sig,"DESeq2_Results_padj_Interaction_Sex-Group_CremTGvsCTR_FvsM.xlsx")

write.table(res_df,"DESeq2_Results_Interaction_Sex-Group_HDAC2_KOvsCTR_FvsM.txt",sep="\t", quote=F,row.names = FALSE)
write.xlsx(res_df,"DESeq2_Results_Interaction_Sex-Group_HDAC2_KOvsCTR_FvsM.xlsx")
write.xlsx(res_df_sig,"DESeq2_Results_padj_Interaction_Sex-Group_HDAC2_KOvsCTR_FvsM.xlsx")










