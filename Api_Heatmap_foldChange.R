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
library(ComplexHeatmap)

#Load data table


  data_foldChange <- read.table('mpare_Heatmap_padj0.05.txt', header=TRUE, row.names=1)
  mat=data_foldChange[,2:4]
  
  mat <- as.matrix(mat)
  

  #Finally a nice explanation https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html 
  Heatmap(mat)
  
  #Define breaks and colours
  library(circlize)
  col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
  col_fun(seq(-3, 3))
  
  #Change breaks with the col function and add title to the legend
  heatmap1<-Heatmap(mat, name = "Log2FC", col = col_fun,
          show_row_names = FALSE,
          cluster_columns = FALSE,
          column_title = "Heatmap of genes affected by Sex:Group interaction \n relative to CTR",
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          column_names_gp= gpar(fontsize=10),
          )
  
  #Save plot
  tiff("Heatmap_of_genes_affected_bySex_Group_interaction_vsCTR.tiff",units = 'in', width = 5, height = 5, res=300)
  draw(heatmap1)
  
  dev.off()


  pheatmap(mat[,2:4],
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           main= "Heatmap of genes affected by Sex:Group interaction",
           border_color = NA,
           scale = 'none',
           #labels_row = mat$GeneSymbol
           show_rownames = FALSE,
           )
  