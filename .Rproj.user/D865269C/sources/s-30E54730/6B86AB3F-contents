


list_LTBi_mirna <- list(miRNA_genes = unique(df_fin$target_symbol), genetodisease_TBi= unique(data_reg_gene2disease_TBi$gene_symbol),enrichment_LTBi = unique(GENENAME_genes_enrichDGN_LTBi$SYMBOL), enrichment_Pulmonary_TB = unique(GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL),enrichDGN_TB_active = unique(GENENAME_genes_enrichDGN_TB_active$SYMBOL))
#list_LTBi_mirna <- list(genetodisease_TBi= unique(data_reg_gene2disease_TBi$gene_symbol),enrichment_LTBi = unique(GENENAME_genes_enrichDGN_LTBi$SYMBOL), enrichment_Pulmonary_TB = unique(GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL),enrichDGN_TB_active = unique(GENENAME_genes_enrichDGN_TB_active$SYMBOL))
list_LTBi_mirna <- list(miRNA_genes = unique(df_fin$target_symbol), 
                        TBi_genes= unique(data_reg_gene2disease_TBi$gene_symbol),
                        LTBi_genes = unique(GENENAME_genes_enrichDGN_LTBi$SYMBOL), 
                        TB_genes = unique(GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL),
                        TB_active_genes = unique(GENENAME_genes_enrichDGN_TB_active$SYMBOL))

## 2-way Venn diagram
Results <- overLapper(list_LTBi_mirna[1:5], type="vennsets")
vennPlot(Results)
library("writexl")
capture.output(Results@vennlist, file = "C:/Users/hp/Desktop/Results.txt")


Results@vennlist[["miRNA_genes_genetodisease_TBi_enrichment_LTBi_enrichment_Pulmonary_TB_enrichDGN_TB_active"]]
Results@vennlist[["enrichment_LTBi"]]
Results@vennlist[["enrichDGN_TB_active"]]
Results@vennlist[["enrichment_LTBi_enrichDGN_TB_active"]]

test_g2 <-df_fin[which( df_fin$target_symbol%in% Results@vennlist[["miRNA_genes_genetodisease_TBi_enrichment_LTBi_enrichment_Pulmonary_TB_enrichDGN_TB_active"]]),]
test_g2

test_g6 <-df_fin[which( df_fin$target_symbol%in% Results@vennlist[["enrichment_LTBi"]]),]
test_g6

test_g9 <-df_fin[which( df_fin$target_symbol%in% Results@vennlist[["enrichDGN_TB_active"]]),]
test_g9

fn <- c(toy2@vennlist[["genetodisease_TBi_enrichment_LTBi_enrichment_Pulmonary_TB"]],datafull_join_getHpa10$SYMBOL,datafull_join_getHpa1d$SYMBOL[1:7])
test_ggg55 <-df_fin[which( df_fin$target_symbol%in% unique(fn)),]
View(test_ggg55)
dim(test_ggg55)
example4@data

mirna.score <-example4@data[which( example4@data$target_symbol%in% test_ggg55$target_symbol),]
mirna.score <-example4@data[which( example4@data$target_symbol%in% c("CSF2" ,    "IL10" ,    "IL1B"  ,   "SLC11A1",  "CD163"  ,  "MAPK1"  ,  "FN1"  ,
                                                                     "CCL2"   ,  "CDH1"   ,  "ANXA1"  ,  "CD209" ,   "TIRAP"  ,  "HLA-DRB1",
                                                                     "HLA-DQA1")),]

mirna.score2 <-  mirna.score[order(mirna.score$score, decreasing = TRUE), ]

#mirna.score <-example4.counts[1:26,][which( example4@data$target_symbol%in% test_ggg55$target_symbol),]
#mirna.score2 <-  mirna.score[order(mirna.score$score, decreasing = TRUE), ]

dim(mirna.score2)
unique(mirna.score2$mature_mirna_id)
mirna.score3 <- mirna.score2[,3]


mirna.score2.counts <- addmargins(table(mirna.score2[, 3:4]))
mirna.score2.counts <- mirna.score2.counts[-nrow(mirna.score2.counts), ]
mirna.score2.counts <- mirna.score2.counts[order(mirna.score2.counts[, 16], decreasing = TRUE), ]

unique(row.names(mirna.score2.counts))



write.csv(mirna.score2.counts, file = 'C:/Users/hp/Desktop/mirna.score2.counts.csv', row.names = FALSE)



library("hpar")
id2 <- c(Results@vennlist[["miRNA_genes_genetodisease_TBi_enrichment_LTBi_enrichment_Pulmonary_TB_enrichDGN_TB_active"]]
,
Results@vennlist[["enrichment_LTBi"]]
,
Results@vennlist[["enrichDGN_TB_active"]]
)
getHpa122 <- getHpa(id2, hpadata = "hpaNormalTissue")
getHpa11223 <-getHpa122[which( getHpa122$Tissue%in% "lung"),]
library("dplyr") 
colnames(getHpa11)[2] <- "SYMBOL"
library(readxl)
tissuesDB <- read_excel("tissuesDB.xlsx")
datafull_join_getHpa11 <- dplyr::full_join(getHpa11[,c("Gene","SYMBOL","Tissue","Cell.type","Level","Reliability")], 
                                           tissuesDB[,c("SYMBOL","TissuesScore_Lung","TissuesScore_ImuneSystem","TissuesScore_Blood")], by='SYMBOL')

library("dplyr")


datafull_join_getHpa10 <- datafull_join_getHpa11 |>
  filter(Level == "High") |>
  arrange(desc("TissuesScore_Lung")) |>
  head()
