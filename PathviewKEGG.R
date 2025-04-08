library(pathview)
library(openxlsx)

#load data as geneID LOG2FC-sample1, Log2FC-sample2 ...
data(gse16873.d)


#load pathway related data
data(demo.paths)

#we can check the full list of pathways 
data(paths.hsa)
head(paths.hsa, 3)

#
i=1
pv.out <- pathview(gene.data = gse16873.d[, 1],
                   pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", 
                   out.suffix = "gse16873", kegg.native = T)

#Can incliude more samples by including the datasets columns. i.e '1:3'
pv.out <- pathview(gene.data = gse16873.d[, 1:3],
                   pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", 
                   out.suffix = "gse16873", kegg.native = T)
#Dataset of species korg
data(korg)
#mus musculus mmu

#----------------Our dataset#####

mypaths <- list("00190", "05012", "05415")

data1 <-read.xlsx("pathviewINPUT_KEGG_mmu00190_OxPh.xlsx", rowNames = TRUE)
data2<-read.xlsx("pathviewINPUT_KEGG_mmu05012_Parkinson.xlsx", rowNames = TRUE)
data3<- read.xlsx("pathviewINPUT_KEGG_mmu05415_DiabCard.xlsx", rowNames = TRUE)

column_range <-c(1,4)

i=3

pv.out <- pathview(gene.data = data3[, column_range],
                   pathway.id = mypaths[i],
                   species = "mmu", 
                   out.suffix = paste(mypaths[i]), kegg.native = T)
