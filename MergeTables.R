library(tidyverse)
library(dplyr)
library(openxlsx)



setwd("W:/01_GenEpi_Projects/Schulte_2025/mergeTables")

df1 <-read.xlsx("DESeq2_Results_Condition_CremTG_F_vs_CremTG_M.xlsx")
df2 <-read.xlsx("DESeq2_Results_Condition_CTR_F_vs_CTR_M.xlsx")
df3 <-read.xlsx("DESeq2_Results_Condition_HDAC2_KO_F_vs_HDAC2_KO_M.xlsx")
df4<-read.xlsx("DESeq2_Results_Condition_TGxKO_F_vs_TGxKO_M.xlsx")
df5<-read.xlsx("DESeq2_Results_Group_F_CremTG_vs_CTR.xlsx")
df6<-read.xlsx("DESeq2_Results_Group_F_HDAC2_KO_vs_CTR.xlsx")
df7<-read.xlsx("DESeq2_Results_Group_F_TGxKO_vs_CremTG.xlsx")
df8<-read.xlsx("DESeq2_Results_Group_F_TGxKO_vs_CTR.xlsx")
df9<-read.xlsx("DESeq2_Results_Group_M_CremTG_vs_CTR.xlsx")
df10<-read.xlsx("DESeq2_Results_Group_M_HDAC2_KO_vs_CTR.xlsx")
df11<-read.xlsx("DESeq2_Results_Group_M_TGxKO_vs_CremTG.xlsx")
df12<-read.xlsx("DESeq2_Results_Group_M_TGxKO_vs_CTR.xlsx")
df13<-read.xlsx("Group_CremTG_vs_CTR_DESeq2_Results.xlsx")
df14<-read.xlsx("Group_HDAC2_KO_vs_CTR_DESeq2_Results.xlsx")
df15<-read.xlsx("Group_TGxKO_vs_CremTG_DESeq2_Results.xlsx")
df16<-read.xlsx("Group_TGxKO_vs_CTR_DESeq2_Results.xlsx")
df17<-read.xlsx("DESeq2_Results_Interaction_Sex-Group_CremTGvsCTR_FvsM.xlsx")
df18<-read.xlsx("DESeq2_Results_Interaction_Sex-Group_TGKOvsCTR_FvsM.xlsx")
df19<-read.xlsx("DESeq2_Results_Interaction_Sex-Group_HDAC2_KOvsCTR_FvsM.xlsx")


library(dplyr)
df_merge<- left_join(df1, df2, by='ENSEMBL') %>%
  left_join(., df3, by='ENSEMBL')  %>%
  left_join(., df4, by='ENSEMBL')  %>%
  left_join(., df5, by='ENSEMBL')  %>%
  left_join(., df6, by='ENSEMBL')  %>%
  left_join(., df7, by='ENSEMBL')  %>%
  left_join(., df8, by='ENSEMBL')  %>%
  left_join(., df9, by='ENSEMBL')  %>%
  left_join(., df10, by='ENSEMBL')  %>%
  left_join(., df11, by='ENSEMBL')  %>%
  left_join(., df12, by='ENSEMBL')  %>%
  left_join(., df13, by='ENSEMBL')  %>%
  left_join(., df14, by='ENSEMBL')  %>%
  left_join(., df15, by='ENSEMBL')  %>%
  left_join(., df16, by='ENSEMBL')  %>%
  left_join(., df17, by='ENSEMBL')  %>%
  left_join(., df18, by='ENSEMBL')  %>%
  left_join(., df19, by='ENSEMBL')  

write.xlsx(df_merge, "20250403_ALL_COMPARISONS_TABLE.xlsx")


colnames(df_merge)[grepl("padj", colnames(df_merge))]

filtered_df_merge <- df_merge %>% filter(!(padj.x >0.1 &
                                             padj.y >0.1 &
                                             padj.x.x >0.1 &
                                             padj.y.y >0.1 &
                                             padj.x.x.x >0.1 &
                                             padj.y.y.y >0.1 &
                                             padj.x.x.x.x >0.1 &
                                             padj.y.y.y.y >0.1 &
                                             padj.x.x.x.x.x >0.1 &
                                             padj.y.y.y.y.y >0.1 &
                                             padj.x.x.x.x.x.x >0.1 &
                                             padj.y.y.y.y.y.y >0.1 &
                                             padj.x.x.x.x.x.x.x >0.1 &
                                             padj.y.y.y.y.y.y.y >0.1 &
                                             padj.x.x.x.x.x.x.x.x >0.1 &
                                             padj.y.y.y.y.y.y.y.y >0.1 &
                                             padj.x.x.x.x.x.x.x.x.x >0.1 &
                                             padj.y.y.y.y.y.y.y.y.y >0.1 &
                                             padj>0.1))

write.xlsx(filtered_df_merge, "20250403_ALL_COMPARISONS_TABLE_padj0.1.xlsx")
                                         

filtered_df_merge_005 <- df_merge %>% filter(!(padj.x >0.05 &
                                                padj.y >0.05 &
                                                padj.x.x >0.05 &
                                                padj.y.y >0.05 &
                                                padj.x.x.x >0.05 &
                                                padj.y.y.y >0.05 &
                                                padj.x.x.x.x >0.05 &
                                                padj.y.y.y.y >0.05 &
                                                padj.x.x.x.x.x >0.05 &
                                                padj.y.y.y.y.y >0.05 &
                                                padj.x.x.x.x.x.x >0.05 &
                                                padj.y.y.y.y.y.y >0.05 &
                                                padj.x.x.x.x.x.x.x >0.05 &
                                                padj.y.y.y.y.y.y.y >0.05 &
                                                padj.x.x.x.x.x.x.x.x >0.05 &
                                                padj.y.y.y.y.y.y.y.y >0.05 &
                                                padj.x.x.x.x.x.x.x.x.x >0.05 &
                                                padj.y.y.y.y.y.y.y.y.y >0.05 &
                                                padj>0.05))
                                                                                      
                  write.xlsx(filtered_df_merge_005, "20250403_ALL_COMPARISONS_TABLE_padj0.05.xlsx")
