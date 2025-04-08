#Get genes for KEGG pathways and plot the heatmap of Log2FoldChanges for these genes for the comparisons of interest
#https://support.bioconductor.org/p/109871/ 


  
#-----------Load libraries


library(org.Mm.eg.db)
library(matrixStats)
library("RColorBrewer")
library(pheatmap)
library(openxlsx)
library(ComplexHeatmap)

#Set working directory
setwd("W:/01_GenEpi_Projects/Schulte_2025/3.Sex&Group/ShinyGo/ComparingGeneSetFoldChanges")
#Load data table
#Select columns for heatmap 

  library(openxlsx)
  data_foldChange <- read.xlsx('FULLTABLE.xlsx')
  #data_heatmap <- data_foldChange[,column_vector]
  data_heatmap <- data_foldChange

#Create Gene names for KEGG pathways
  #https://support.bioconductor.org/p/109871/
  
  library("KEGGREST")
  library(org.Mm.eg.db)
  library(tidyverse)
  browseVignettes("KEGGREST")
  
  #Get KEGG pathways and Entrex gene ids
  #transform to a tibble data.frame with tydiverse::tibble()
  mmu_path_eg <- keggLink("pathway", "mmu") %>%
    tibble(pathway =., eg=sub("mmu:", "", names(.)))
    
  #annotate with SYMBOL and ENSEMBL identifiers
  mmu_kegg_anno <- mmu_path_eg %>%
    mutate( GeneSymbol = mapIds(org.Mm.eg.db, eg, "SYMBOL", "ENTREZID"),
            ensembl = mapIds(org.Mm.eg.db, eg, "ENSEMBL", "ENTREZID"))
  mmu_kegg_anno$pathway <-str_remove(mmu_kegg_anno$pathway, "path:")
  
  #Create pathway description from KEGG 
  mmu_pathways <- keggList("pathway", "mmu") %>%
    tibble (pathway = names(.), description =.)
  
#Extract list of genes for specific pathway
  my_pathways = c("mmu05012", #Parkinson's disease
                  "mmu05415", #Diabetic cardiomyopathy
                  "mmu00190") #Oxidative phosphorylation
  
  #Extract list of genes for each pathway
  
#####KEGG Parkinson disease
  
  mmu_ParkDis <- mmu_kegg_anno [mmu_kegg_anno$pathway == "mmu05012",]
  mmu_ParkDis_mydata <- left_join( mmu_ParkDis, data_heatmap, by ="GeneSymbol" )
  #write table to keep
  #write.xlsx(mmu_ParkDis_mydata, "KEGG_mmu05012_Parkinson.xlsx")
  
  #Create table for heatmap
  #Drop missing values for ENSEMBL - genes that are missing from mydata 
  #tidyr::drop_na()
  mmu_ParkDis_mydata_hmap <- mmu_ParkDis_mydata %>% drop_na(ENSEMBL)
  
  #Keep only genes where one of the comparisons is statistically significant 
  mmu_ParkDis_mydata_hmap <- subset( mmu_ParkDis_mydata_hmap, mmu_ParkDis_mydata_hmap$CremTG_padj< 0.05|
                                       mmu_ParkDis_mydata_hmap$TGKO_padj <0.05 |
                                       mmu_ParkDis_mydata_hmap$TGKOvsCremTG_padj <0.05 )
  
  #for subsetting 
  column_vector <- c(3,6,9)
  mmu_ParkDis_mydata_hmap <- mmu_ParkDis_mydata_hmap[,column_vector]
  
  
 
  
######KEGG Diabetic cardiomyopathy
  
  mmu_DiabCard <- mmu_kegg_anno [mmu_kegg_anno$pathway == "mmu05415",]
  mmu_DiabCard_mydata <- left_join( mmu_DiabCard, data_heatmap, by ="GeneSymbol" )
  #write table to keep
 # write.xlsx(mmu_DiabCard_mydata, "KEGG_mmu05415_DiabCard.xlsx")
 
  
  #Drop missing values for ENSEMBL - genes that are missing from mydata 
  #tidyr::drop_na()
  mmu_DiabCard_mydata_hmap <- mmu_DiabCard_mydata %>% drop_na(ENSEMBL)
  
  #Keep only genes where one of the comparisons is statistically significant 
  mmu_DiabCard_mydata_hmap <- subset( mmu_DiabCard_mydata_hmap, mmu_DiabCard_mydata_hmap$CremTG_padj< 0.05|
                                       mmu_DiabCard_mydata_hmap$TGKO_padj <0.05 |
                                       mmu_DiabCard_mydata_hmap$TGKOvsCremTG_padj <0.05 )
  
  #for subsetting 
  column_vector <- c(3,6,9)
  mmu_DiabCard_mydata_hmap <- mmu_DiabCard_mydata_hmap[,column_vector]
  
  
######KEGG  Oxidative Phosphorylation
  
  mmu_OxPh <- mmu_kegg_anno [mmu_kegg_anno$pathway == "mmu00190",]
  mmu_OxPh_mydata <- left_join( mmu_OxPh, data_heatmap, by ="GeneSymbol" )
  
  #write table to keep
  #write.xlsx(mmu_OxPh_mydata, "KEGG_mmu00190_OxPh.xlsx")
  
  
  #Drop missing values for ENSEMBL - genes that are missing from mydata 
  #tidyr::drop_na()
  mmu_OxPh_mydata_hmap <- mmu_OxPh_mydata %>% drop_na(ENSEMBL)
  
 
  #Keep only genes where one of the comparisons is statistically significant 
  mmu_OxPh_mydata_hmap <- subset( mmu_OxPh_mydata_hmap, mmu_OxPh_mydata_hmap$CremTG_padj< 0.05|
                                        mmu_OxPh_mydata_hmap$TGKO_padj <0.05 |
                                        mmu_OxPh_mydata_hmap$TGKOvsCremTG_padj <0.05 )
  
  #for subsetting 
  column_vector <- c(3,6,9)
  mmu_OxPh_mydata_hmap <- mmu_OxPh_mydata_hmap[,column_vector]
  
  
  
  
    
    
    
    
    
  
  
  
  
  
#########HEATMAPS###################### 
  
  { 
  mat <- as.matrix(mmu_ParkDis_mydata_hmap[,2:3])
  
  #Finally a nice explanation https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html 
  Heatmap(mat)
  
  #Define breaks and colours
  #library(circlize)
  #col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  #col_fun(seq(-3, 3))
  
  #Change breaks with the col function and add title to the legend
  heatmap1<-Heatmap(mat, name = "Log2FC", 
                    #col = col_fun,
          show_row_names = FALSE,
          cluster_columns = FALSE,
          column_title = paste("Heatmap of genes in KEGG pathway 'Parkinson's disease' \n relative to CTR"),
          column_title_gp = gpar(fontsize = 10, fontface = "bold"),
          column_names_gp= gpar(fontsize=8),
          )
  
  #Save plot
  tiff("Heatmap_of_genes_KEGGPark_vsCTR.tiff",units = 'in', width = 5, height = 5, res=300)
  draw(heatmap1)
  
  dev.off()
  }
  
  { 
    mat <- as.matrix(mmu_DiabCard_mydata_hmap[,2:3])
    
    #Finally a nice explanation https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html 
    Heatmap(mat)
    
    #Define breaks and colours
    #library(circlize)
    #col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
    #col_fun(seq(-3, 3))
    
    #Change breaks with the col function and add title to the legend
    heatmap1<-Heatmap(mat, name = "Log2FC", 
                      #col = col_fun,
                      show_row_names = FALSE,
                      cluster_columns = FALSE,
                      column_title = paste("Heatmap of genes in KEGG pathway 'Diabetic Cardiomyopathy' \n relative to CTR"),
                      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      column_names_gp= gpar(fontsize=8),
    )
    
    #Save plot
    tiff("Heatmap_of_genes_KEGGDiabCard_vsCTR.tiff",units = 'in', width = 5, height = 5, res=300)
    draw(heatmap1)
    
    dev.off()
  }

  { 
    mat <- as.matrix(mmu_OxPh_mydata_hmap[,2:3])
    
    #Finally a nice explanation https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html 
    Heatmap(mat)
    
    #Define breaks and colours
    #library(circlize)
    #col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
    #col_fun(seq(-3, 3))
    
    #Change breaks with the col function and add title to the legend
    heatmap1<-Heatmap(mat, name = "Log2FC", 
                      #col = col_fun,
                      show_row_names = FALSE,
                      cluster_columns = FALSE,
                      column_title = paste("Heatmap of genes in KEGG pathway 'Oxidative Phosphorylation \n relative to CTR"),
                      column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                      column_names_gp= gpar(fontsize=8),
    )
    
    #Save plot
    tiff("Heatmap_of_genes_KEGGOxPh_vsCTR.tiff",units = 'in', width = 5, height = 5, res=300)
    draw(heatmap1)
    
    dev.off()
  }
