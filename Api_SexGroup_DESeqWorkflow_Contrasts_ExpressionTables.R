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

#Comprehensive directory
Mum_dir='W:/01_GenEpi_Projects/Schulte_2025/3.Sex&Group/'
#Create metrics data.frame to append result metrics of each comparison to
Summary_table <- data.frame(matrix(ncol=4, nrow=0))
Summary_table_cols <- c("Comparison", "padj_sign_genes","Upregulated", "Downregulated")
colnames(Summary_table) <- Summary_table_cols


#-----------Define Variables
##Working directory
  wdir= Mum_dir
##Control group to use as reference level for the DESeq Object
CTRL_level='CTR'

## FDR cutoff for results() and results summary 
  #padj threshold for statistical significance (for alpha=)
  alpha_var =0.05
  #Log2FoldChange threshold for the Wald test (for lfcThreshold=). Default altHypothesis tests for
  #LFC > or < of absolute value specified by lfcThreshold=
  lfcThreshold_var=0
  


## Volcano plot variables
  #yaxis upper limit
  #ylim= 10

  #xaxis lower limit
  #left_xlim=-8
  #xaxis upper limit
  #right_xlim=8

  #cutoffs for significance
  #Fold change threshold in absolute value
  FC_thr=0
  #pvalue thresholds
  pval_thr=0.05

##


#-----------Load data+ check data format

  #Set working directory
  setwd(paste(wdir))

  #Input data: count matrix generated with Htseq-count where the column names correspond to sample names, in the same order
   #as the row names in the cov matrix
  count_matrix <- read.table('CountMatrix.txt', header=TRUE, row.names=1)

  #Matrix containing the metadata for the samples, row names correspond to column names in count_matrix
  #Contains the 'Group' variable to be used in the experimental design and that includes the reference level for comparisons
  cov <- read.table('Groups.txt', header=TRUE, row.names=1)

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
                            design= ~ Sex + Group)
dds



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
  dds$Group <-relevel (dds$Group, ref = paste(CTRL_level))
  dds$Sex <-relevel (dds$Sex, ref="M")

#------------Run DESeq function
  #In this step the data is normalised by the size factor
  dds <-DESeq(dds)
  resultsNames(dds)
  
  res <- results(dds)
  res

  #Explore results
  summary(res)

  
###################################Save normalised and transformed datasets for this specific analysis design
  ########Annotate count tables with GENEID and saved normalised data, transformed data. 
  
########--------write table of normalised counts
  norm_counts<-counts(dds, normalized=TRUE)
  norm_counts<-as.data.frame(norm_counts)
  ann_norm_counts <- merge (norm_counts, GENEID, by=0)
  
  #reorder columns 
  #library(dplyr)
  #remove redundant ENSEMBLIDs column
  ann_norm_counts <-ann_norm_counts[,-1]
  
  #Reorder columns to have ENSEMBLIDs, GeneSymbol and then samples 
  ann_norm_counts <- ann_norm_counts %>%
    dplyr::select(ENSEMBL, GeneSymbol, everything())
  
  #write data to a txt file 
  write.table(ann_norm_counts,'normalised_counts.txt',sep="/t", quote=F,row.names = FALSE)
  write.xlsx(ann_norm_counts, 'normalised_counts.xlsx')
  
  
########  save vst data
  vsd <- vst(dds, blind=TRUE)
  vsd <- assay(vsd)
  vsd_df <- as.data.frame(vsd)
  vsd_df <- merge (vsd_df, GENEID, by=0)
  
  #reorder columns 
  #library(dplyr)
  #remove redundant ENSEMBLIDs column
  vsd_df <-vsd_df[,-1]
  
  #Reorder columns to have ENSEMBLIDs, GeneSymbol and then samples 
  vsd_df <- vsd_df %>%
    dplyr::select(ENSEMBL, GeneSymbol, everything())
  
  #write data to a txt file 
  write.table(vsd_df,'vst_counts.txt',sep="/t", quote=F,row.names = FALSE)
  write.xlsx(vsd_df, 'vst_counts.xlsx')
  

  
####CHeck that the two transformations have worked as expected
#From the vignette: 'The regularized log to a lesser extent, while for the variance stabilized data
    #the standard deviation is roughly constant along the whole dynamic range.
    #Note that the vertical axis in such plots is the square root of the variance over all samples, 
    #so including the variance due to the experimental conditions. 
    #While a flat curve of the square root of variance over the mean may seem like the goal of such 
    #transformations, this may be unreasonable in the case of datasets with many true differences due 
    #to the experimental conditions.
  
  library(vsn)
  tiff("vsd_stdevOverMean.tiff",units = 'in', width = 5, height = 5, res=300)
  plot_vsd <- meanSdPlot(vsd)
  dev.off()
 
################################################################################################
############# Differential Expression Analysis - contrasts, volcano plots ... ##################

  #Comparisons loop
  #Levels variable to specify the level of the comparison for the denominator
  Levels <- list(c("CTR"),
                 c("CTR"),
                 c("CTR"),
                 c("CremTG"))
  
  #Create a list with character vectors for each comparison to be assigned to the 'comparison' term in 'result' function
    # c('VariableName', 'NumeratorTermForFoldChange', 'DenominatorTermForFoldChange')
  Comparisons <- list(c("Group","HDAC2_KO","CTR"),
                      c("Group", "CremTG", "CTR"),
                      c("Group", "TGxKO", "CTR"),
                      c("Group","TGxKO","CremTG"))

  number_comparisons <-length(Comparisons)
  #Select element from 'Comparisons'
  #i=4
  
  for (i in 1:number_comparisons) 
    {
    
    #set the factor level
    dds$Group <-relevel (dds$Group, ref = Levels[[i]])
    
   
    #We need to rerun the DESeq function only if we change the level
    if (i==4) 
    {
     #------------Run DESeq function
    #In this step the data is normalised by the size factor
    dds <-DESeq(dds)
    resultsNames(dds)
    }
    
    #Calculate fold change and p-value with Wald test
    #By specifying lfcThreshold we will test for genes that have a LFC greater of the absolute value specified
      #It is preferable to set the LFC threshold during this step, rather then in downstream filtering, since this way 
      #the actual number of tests is captured for correcting for multiple testing with the BH method. 
    res<-results(dds, contrast= Comparisons[[i]],
                 lfcThreshold = lfcThreshold_var) 
    
    #Create label for the comparison
    term1 <- Comparisons[[i]][[1]]
    term2 <- Comparisons[[i]][[2]]
    term3 <- Comparisons[[i]][[3]]
    comparison_lab <-paste(term1,'_',term2,'_vs_', term3, sep='')
    #Create directory of results for the comparison and assign it to out_dir  
    #dir.create(comparison_lab)
    #out_dir <- paste0(comparison_lab,"/")
    
    ####Save summary of results to txt 
    sink(paste0(comparison_lab,"Summary_DESeq2_Results.txt"))
    summary(res, alpha = alpha_var)
    sink(NULL)
    
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
    
    
    #Calculate the shrunken Fold changes as calculated by the DESeq function
    #The default method used by DESeq is apeglm
    #These shrunken estimates are used for input in the Wald test for statistical significance
    #The authors show how these shrunken LFC are more accurate for subsetting significant hits (improving reproducibility)
    #The coefficient is the variable of interest for which we want to calculate the shrunk fold changes 
      #as present in resultsNames(dds)
    library(apeglm)
    resLFCSh <-lfcShrink(dds, coef=comparison_lab, type ='apeglm', lfcThreshold = lfcThreshold_var)
    resLFCSh_df <- as.data.frame(resLFCSh)
    #change Column names 
    colnames(resLFCSh_df)[1]<-"Shr_baseMean"
    colnames(resLFCSh_df)[2]<-"Shr_LFC"
    colnames(resLFCSh_df)[3]<-"Shr_lfcSE"
    colnames(resLFCSh_df)[4]<-"Shr_svalue"
    #Create ENSEMBL column from rownames
    resLFCSh_df$ENSEMBL <- rownames(resLFCSh_df)
    head(resLFCSh_df)
    #Create subset for annotation
    cols=c(2,6)
    resLFCSh_df_sub <-resLFCSh_df[,cols]
    head(resLFCSh_df_sub)

    
####FOR HEATMAP: Create a subset dataset with only genes with padj <pval_thr sorted by decreasing |Shrunk_LFC|
    #Keep only genes with padj <pval_thr
    res_df_sig<- res_df[res_df$padj <= pval_thr,]
    #Remove rows where padj = missing 
    res_df_sig <- res_df_sig %>% drop_na(padj)

    #Arrange by decreasing absolute value of LFC
    res_df_sig <- res_df_sig %>% arrange(desc(abs(log2FoldChange)))
    #Add ShrunkLFC column
    res_df_sig_shrLFC<- left_join(res_df_sig, resLFCSh_df_sub, by='ENSEMBL')
    head(res_df_sig_shrLFC)
    #Arrange by decreasing absolute value of Shrunk_LFC
    res_df_sig_shrLFC <- res_df_sig_shrLFC %>% arrange(desc(abs(Shr_LFC)))
    head(res_df_sig_shrLFC)
    #Create Vector with TOP50 genes 
    Top50shrLFC <- res_df_sig_shrLFC$ENSEMBL[1:50]
    
    
    
    
    ###plot-MA 
    plotMA<- plotMA (res, ylim=c(-2,2))
    tiff(paste0(comparison_lab,"_plotMA.tiff"),units = 'in', width = 5, height = 5, res=300)
    print(plotMA)
    dev.off()
    
    plotMA_shrunk <-plotMA(resLFCSh,ylim=c(-2,2))
    tiff(paste0(comparison_lab,"_plotMA_shrunk.tiff"),units = 'in', width = 5, height = 5, res=300)
    print(plotMA_shrunk)
    dev.off()
    
  ####Prepare table for saving results
    #Add shrunkLFC to result table 
    res_df<-left_join(res_df, resLFCSh_df_sub, by='ENSEMBL')
    #Sort data by increasing adj pvalues    
    library(tidyverse)
    res_df <-  res_df %>%
      arrange(padj)  # arrange in ascending order, NAs are put at the bottom of the file
    head(res_df)
    #Create a vector with these top significant genes (to be used to subset dataset for heatmap)
    #Use ENSEMBLID because they are unique to avoid issues
    head(res_df)
    Top50padj <- res_df$ENSEMBL[1:50]
    
    
    write.table(res_df,paste0(comparison_lab,"_DESeq2_Results.txt"),sep="/t", quote=F,row.names = FALSE)
    write.xlsx(res_df,paste0(comparison_lab,"_DESeq2_Results.xlsx"))
    write.xlsx(res_df_sig_shrLFC,paste0(comparison_lab,"_DESeq2_Results_padj",pval_thr,".xlsx"))
    
    #Save metrics for summary table
    num_sign_genes= nrow(res_df_sig_shrLFC)
    UPs<- res_df_sig_shrLFC[res_df_sig_shrLFC$log2FoldChange>0,]
    DOWNs<- res_df_sig_shrLFC[res_df_sig_shrLFC$log2FoldChange<0,]
    Upreg_genes= nrow(UPs)
    Downreg_genes=nrow(DOWNs)
    
    #Add metrics to summary table
    #Summary_table <- data.frame(matrix(ncol=4, nrow=0))
    #Summary_table_cols <- c("Comparison", "padj_sign_genes","Upregulated", "Downregulated")
   info_list <- list(comparison_lab, num_sign_genes, Upreg_genes, Downreg_genes)
    Summary_table[nrow(Summary_table)+1,] <- info_list
    head(Summary_table)
    
########  Volcano plot
    #Loading relevant libraries
    library (tidyverse)
    library(RColorBrewer)
    library(ggrepel)
    

    #assign dataset 
    df = res_df
    comparison = comparison_lab
    Title <- paste('Volcano plot for',Comparisons[[i]][[2]],"vs",Comparisons[[i]][[3]])
    
    # Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
    df$diffexpressed <- "NO"
    # if log2Foldchange > FC_thr and pvalue < 0.05, set as "UP"
    df$diffexpressed[df$log2FoldChange > FC_thr & df$padj < pval_thr] <- "UP"
    # if log2Foldchange < FC_thr and pvalue < 0.05, set as "DOWN"
    df$diffexpressed[df$log2FoldChange < -FC_thr & df$padj < pval_thr] <- "DOWN"
    # Explore a bit
    head(df[order(df$padj) & df$diffexpressed == 'DOWN', ])
    
    # Create a new column "delabel" that will contain the name of the top 50 differentially expressed genes (NA in case they are not)
    df$delabel <- ifelse(df$GeneSymbol %in% head(df[order(df$padj), "GeneSymbol"], 50), df$GeneSymbol, NA)
  
    #Add specific genes of interest
    #df$delabel[df$GENEID=='VASH1']<-'VASH1'
    #df$delabel[df$GENEID=='GUCY1A1']<-'GUCY1A1'
    
    myvolcanoplot <- ggplot(data = df, aes(x = log2FoldChange,
                                           y = -log10(pvalue),
                                           col = diffexpressed,
                                           label = delabel)) +
      geom_vline(xintercept = c(0),
                 col = "gray", 
                 linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), 
                 col = "gray", 
                 linetype = 'dashed') +
      geom_point(size = 0.3) +
      scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3"), # to set the colours of our variable
                         labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
      #coord_cartesian(ylim = c(0, ylim), xlim = c(left_xlim, right_xlim)) + # since some genes can have minuslog10padj of inf, we set these limits
      labs(color = '', #legend_title,<br />
           x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) +
      scale_x_continuous(breaks = seq(-100, 100, 2)) + # to customise the breaks in the x axis
      ggtitle(paste(Title)) + # Plot title
      
      geom_text_repel(max.overlaps = Inf,
                      size=4,
                      segment.size=0.2)  
    # Open the file that will contain your plot (the name is up to you)
    #pdf(file = paste('VolcanoPlot_',comparison, '.pdf', sep=''), width = 8, height = 12) # you can change the size of the output file
    #pdf(file = paste('VolcanoPlot_',comparison, '.pdf', sep=''), width = 8, height = 12) # you can change the size of the output file
    tiff(paste0(comparison_lab,'_VolcanoPlot.tiff', sep=''),units = 'in', width = 8, height = 12, res=300)
    
    # Execute the plot
    print(myvolcanoplot)
    # Close the file that will contain the plot
    dev.off()
    
    
########  Heatmap with vst data - Top50 genes by padj
    
    #Create subset of transformed dataset for the samples of interest 
    library(dplyr)
    #Identify samples that match the groups of the contrast
    cov_subset <- cov %>% filter(Group == paste(Comparisons[[i]][[2]]) | Group == paste(Comparisons[[i]][[3]]))
    #Create subset of transformed dataset only for those samples
    #Create vector for subsetting the dataset for the samples of interest
    Subset_v <- row.names(cov_subset)
    #Add the names of 'ENSEMBL' and 'GeneSymbol' columns to the vector for subsetting 
    #ENSEMBL in position 1
    #GeneSymbol in position 2
    Subset_v <-append ('GeneSymbol', Subset_v)
    Subset_v <-append ('ENSEMBL', Subset_v)
    Subset_v
    
    #Subset dataset of transformed data for heatmap
    vsd_df_subset<- subset(vsd_df, select = Subset_v)
    head(vsd_df_subset)
    #Select only the data for the Top50 most significant genes sorted by padj (Previously created vectr 'Top50padj')
    vsd_df_subset_Top50 <- vsd_df_subset[vsd_df_subset$ENSEMBL %in% Top50padj,]
    
    #Replace missing GeneSymbol values for original ENSEMBL IDs 
    library(dplyr)
    vsd_df_subset_Top50<-vsd_df_subset_Top50 %>%
      mutate(GeneSymbol = coalesce(GeneSymbol, ENSEMBL)) # works
    
    #Reorder the dataset to match the order of Top50adj
    vsd_df_subset_Top50<-vsd_df_subset_Top50[match(Top50padj, vsd_df_subset_Top50$ENSEMBL),]
    head(vsd_df_subset_Top50[,1:2])
    head(Top50padj)
    
    #Generate dataframe for annotation of the heatmap with the 'Group' variable
    df_Group <- as.data.frame(colData(dds)[,c("Group")])
    colnames(df_Group) <- c('Group')
    rownames(df_Group)<-rownames(colData(dds))
    

    
    #Plot heatmap samples clustered for the first 50 most significant genes sorted by padj (the same genes labelled in the volcano plot)
    vst_heatmap_TOP<- pheatmap(vsd_df_subset_Top50[,3:ncol(vsd_df_subset_Top50)], 
             cluster_rows=FALSE, 
             show_rownames=TRUE,
             labels_row = vsd_df_subset_Top50$GeneSymbol,# to use GeneSymbol to label the rows in the heatmap
             cluster_cols=TRUE, 
             annotation_col=df_Group, 
             annotation_legend=TRUE, 
             #clustering_method ="ward.D",
             fontsize = 5.5,
             main= "Heatmap of TOP50 genes by adjusted-pvalue")
    #Save plot
    tiff(paste0(comparison_lab,"_Heatmap_TOP50genes_bypadj_vstdata.tiff"),units = 'in', width = 5, height = 5, res=300)
    print(vst_heatmap_TOP)
    dev.off()
  
########  Heatmap with vst data - Genese sorted by shrunkFoldChange values
    
    #Create subset of transformed dataset for the samples of interest 
    library(dplyr)
    #Identify samples that match the groups of the contrast
    cov_subset <- cov %>% filter(Group == paste(Comparisons[[i]][[2]]) | Group == paste(Comparisons[[i]][[3]]))
    #Create subset of transformed dataset only for those samples
    #Create vector for subsetting the dataset for the samples of interest
    Subset_v <- row.names(cov_subset)
    #Add the names of 'ENSEMBL' and 'GeneSymbol' columns to the vector for subsetting 
    #ENSEMBL in position 1
    #GeneSymbol in position 2
    Subset_v <-append ('GeneSymbol', Subset_v)
    Subset_v <-append ('ENSEMBL', Subset_v)
    Subset_v
    
    #Subset dataset of transformed data for heatmap
    vsd_df_subset<- subset(vsd_df, select = Subset_v)
    head(vsd_df_subset)
    
    #Select only the data for the Top50 significant genes sorted by highest |shrunkLFC|
    #Create the vector with the names of these genes
    #add shrunkLFC column to the result dataset 
    vsd_df_subset_Top50shrLFC <- vsd_df_subset[vsd_df_subset$ENSEMBL %in% Top50shrLFC,]
    #Replace missing GeneSymbol values for original ENSEMBL IDs 
    library(dplyr)
    vsd_df_subset_Top50shrLFC<-vsd_df_subset_Top50shrLFC %>%
      mutate(GeneSymbol = coalesce(GeneSymbol, ENSEMBL)) # works
    
    #Generate dataframe for annotation of the heatmap with the 'Group' variable
    df_Group <- as.data.frame(colData(dds)[,c("Group")])
    colnames(df_Group) <- c('Group')
    rownames(df_Group)<-rownames(colData(dds))
    
    #Reorder the dataset to match the order of Top50shrLFC
    vsd_df_subset_Top50shrLFC<-vsd_df_subset_Top50shrLFC[match(Top50shrLFC, vsd_df_subset_Top50shrLFC$ENSEMBL),]
    head(vsd_df_subset_Top50shrLFC[,1:2])
    head(Top50shrLFC)
    
    #Plot heatmap samples clustered for the first 50 most significant genes sorted by shrLFC (the same genes labelled in the volcano plot)
    vst_heatmap_TOP_shrLFC<- pheatmap(vsd_df_subset_Top50shrLFC[,3:ncol(vsd_df_subset_Top50shrLFC)], 
                               cluster_rows=FALSE, 
                               show_rownames=TRUE,
                               labels_row = vsd_df_subset_Top50shrLFC$GeneSymbol,# to use GeneSymbol to label the rows in the heatmap
                               cluster_cols=TRUE, 
                               annotation_col=df_Group, 
                               annotation_legend=TRUE, 
                               #clustering_method ="ward.D",
                               fontsize = 5,
                               main= "Heatmap of TOP50 genes by shrunkLFC")
    #Save plot
    tiff(paste0(comparison_lab,"_Heatmap_TOP50genes_byshrLFC_vstdata.tiff"),units = 'in', width = 5, height = 5, res=300)
    print(vst_heatmap_TOP_shrLFC)
    dev.off()
    
    
    
    #end
    }

  #Write summary table
  
  write.xlsx(Summary_table, paste0(Mum_dir,'SummaryTable.xlsx'))
    
    




