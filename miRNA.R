## ---- include=FALSE, echo=FALSE-----------------------------------------------
# date: "`r doc_date()`"
# "`r pkg_ver('BiocStyle')`"
# <style>
#     pre {
#     white-space: pre !important;
#     overflow-y: scroll !important;
#     height: 50vh !important;
#     }
# </style>

## ----annotate, echo=FALSE-------------------------------------------------------------------------
library(knitr)
options(width=100)
opts_chunk$set(echo       = TRUE,
               message    = TRUE,
               warning    = TRUE,
               eval       = TRUE)

## ----multimir_dbInfoVersions----------------------------------------------------------------------
library(multiMiR)
db.ver = multimir_dbInfoVersions()
db.ver

## ----multimir_switchDBVersion, echo=TRUE----------------------------------------------------------
vers_table <- multimir_dbInfoVersions()
vers_table

multimir_switchDBVersion(db_version = "2.0.0")

curr_vers  <- vers_table[1, "VERSION"]  # current version
multimir_switchDBVersion(db_version = curr_vers)

## ----multimir_dbTables----------------------------------------------------------------------------
db.tables = multimir_dbTables()
db.tables

## ----multimir_dbInfo------------------------------------------------------------------------------
db.info = multimir_dbInfo()
db.info

## ----multimir-tabletype---------------------------------------------------------------------------
predicted_tables()
validated_tables()
diseasedrug_tables()
reverse_table_lookup("targetscan")

## ----multimir_dbCount-----------------------------------------------------------------------------
db.count = multimir_dbCount()
db.count
apply(db.count[,-1], 2, sum)

## ----list_multimir--------------------------------------------------------------------------------
miRNAs   = list_multimir("mirna", limit = 10)
genes    = list_multimir("gene", limit = 10)
drugs    = list_multimir("drug", limit = 10)
diseases = list_multimir("disease", limit = 10)
# executes 2 separate queries, giving 20 results
head(miRNAs)
head(genes)
head(drugs)
head(diseases)

## ---- biocworkflow, eval=TRUE---------------------------------------------------------------------
library(edgeR)
library(multiMiR)


## ----Example1-------------------------------------------------------------------------------------
# The default is to search validated interactions in human
example1 <- get_multimir(mirna = c('hsa-miR-4528',"hsa-miR-6844"), summary = TRUE)
names(example1)
# Check which types of associations were returned
table(example1@data$type)
# Detailed information of the validated miRNA-target interaction
head(example1@data)
# Which interactions are supported by Luciferase assay?
example1@data[grep("Luciferase", example1@data[, "experiment"]), ]
example1@summary[example1@summary[,"target_symbol"] == "KRAS",]




## ----Example4_part1-------------------------------------------------------------------------------
example4 <- get_multimir(org     = 'hsa',
                         target  = unique(results_gene2dis_data$gene_symbol),
                         table   = 'predicted',
                         summary = TRUE,
                         predicted.cutoff.type = 'n',
                         predicted.cutoff      = 500000)
View(example4@data)
## ----Example4_part2-------------------------------------------------------------------------------
example4.counts <- addmargins(table(example4@summary[, 2:3]))
example4.counts <- example4.counts[-nrow(example4.counts), ]
example4.counts <- example4.counts[order(example4.counts[, 237], decreasing = TRUE), ]
#write.csv(example4.counts, file = 'C:/Users/hp/Desktop/example4.counts.csv', row.names = FALSE)

head(example4.counts)
View(example4.counts)

d <- example4.counts
names <- rownames(example4.counts)
rownames(d) <- NULL
data <- cbind(names,d)
data

View(example4.counts[1:26,])
data_50 <- data[1:26,]
#data_505 <- data_50[c("names","Sum"),]
#View(data_505)


anno_subset <- example4@data[which( example4@data$mature_mirna_id%in% row.names(example4.counts[1:26,])),]
anno_subset <- anno_subset[,3:4]
length(unique(anno_subset$mature_mirna_id))
length(unique(anno_subset$target_symbol))
dim(anno_subset)
library(dplyr)
df_fin <- anno_subset %>% distinct()
df_fin
dim(df_fin)
View(df_fin)
length(unique(df_fin$mature_mirna_id))
length(unique(df_fin$target_symbol))
unique(df_fin$mature_mirna_id)
list_LTBi_mirna <- list(miRNA_genes = unique(df_fin$target_symbol), genetodisease_TBi= unique(data_reg_gene2disease_TBi$gene_symbol),enrichment_LTBi = unique(GENENAME_genes_enrichDGN_LTBi$SYMBOL), enrichment_Pulmonary_TB = unique(GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL),enrichDGN_TB_active = unique(GENENAME_genes_enrichDGN_TB_active$SYMBOL))
#list_LTBi_mirna <- list(genetodisease_TBi= unique(data_reg_gene2disease_TBi$gene_symbol),enrichment_LTBi = unique(GENENAME_genes_enrichDGN_LTBi$SYMBOL), enrichment_Pulmonary_TB = unique(GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL),enrichDGN_TB_active = unique(GENENAME_genes_enrichDGN_TB_active$SYMBOL))

## 2-way Venn diagram
Results <- overLapper(list_LTBi_mirna[1:5], type="vennsets")
vennPlot(Results)
Results@vennlist[["miRNA_genes_genetodisease_TBi_enrichment_LTBi_enrichment_Pulmonary_TB_enrichDGN_TB_active"]]
Results@vennlist[["enrichment_LTBi"]]
Results@vennlist[["enrichDGN_TB_active"]]
Results@vennlist[["genetodisease_TBi"]]
# install.packages('venneuler')

p1 <- overLapper(list_LTBi[1:2], type="vennsets")
vennPlot(p1)
p2 <- overLapper(list_LTBi[c(1,3)], type="vennsets")
vennPlot(p2)
geneltbis <- p2@vennlist[["miRNA_genes_enrichment_LTBi"]]
unique(geneltbis)

p3 <- overLapper(list_LTBi[c(1,4)], type="vennsets")
vennPlot(p3)
p4 <- overLapper(list_LTBi[c(1,5)], type="vennsets")
vennPlot(p4)
#############################################

list_LTBi_mirna_set2 <-df_fin[which( df_fin$target_symbol%in% geneltbis),]
length(unique(list_LTBi_mirna_set2$target_symbol))
length(unique(list_LTBi_mirna_set2$mature_mirna_id))


test_g <-df_fin[which( df_fin$mature_mirna_id%in% "hsa-miR-424-5p"),]
test_g
test_ggg <-df_fin[which( test_g$target_symbol%in% geneltbis),]
unique(test_ggg$target_symbol)
length(unique(test_ggg$target_symbol))

datadddd <- gene2disease( gene = geneltbis[1:10], verbose = TRUE)
results_datadddd <- datadddd(data3)
View(results_datadddd)
###############################################
library(VennDiagram)
library(purrr)
library(RVenn)
library(ggplot2)
toy = Venn(list_LTBi)
overlap(toy)
#####Draw the Venn Diagram
ggvenn(toy, slice = c(1, 2,3))









library("rbioapi")
mirs <- row.names(example4.counts[1:26,])
mieaa_all <- rba_mieaa_enrich(test_set = mirs,
                              mirna_type = "mature",
                              test_type = "ORA",
                              species = 9606)
mieaa_all


data2 <- example4.counts                                           # Duplicate example data
data2 <- tibble::rownames_to_column(data2, "row_names") # Apply rownames_to_column
data2  

data3 <- example4.counts                                           # Duplicate example data
data3 <- setDT(data3, keep.rownames = TRUE)[]           # Apply setDT function
data3    

## ----annodbi--------------------------------------------------------------------------------------
# On example4's result
columns(example4)
head(keys(example4))
keytypes(example4)
mykeys <- keys(example4)[1:479]
head(select(example4, keys = mykeys, 
            columns = c("database","mature_mirna_acc","mature_mirna_id" ,"score","target_entrez","target_symbol","type")))

View(select(example4, keys = mykeys, 
            columns = c("database","mature_mirna_acc","mature_mirna_id" ,"score","target_entrez","target_symbol","type")))

ff <- select(example4, keys = mykeys, 
       columns = c("database","mature_mirna_acc","mature_mirna_id" ,"score","target_entrez","target_symbol"))
tr <- ff[ff[,"score"] >80.00,]
View(tr)
unique(tr$target_symbol)
# Search by gene on example4's result
columns(example4)
keytypes(example4)
head(keys(example4, keytype = "target_entrez"))
mykeys <- keys(example4, keytype = "target_entrez")[1]
head(select(example4, keys = mykeys, keytype = "target_entrez",
            columns = c("database", "target_entrez", "score")))



library(MiRSEA)
library(miRNApath)






toydd <-AnnotationDbi::select(hgu133plus2.db,keys = unique(list_LTBi_mirna_set2$target_symbol),columns = c("SYMBOL", "ENTREZID","ENSEMBL"),keytype = "SYMBOL")

library("hpar")

## ----getHpa-------------------------------------------------------------------
id <- toydd$ENSEMBL
#Normal tissue data: Expression profiles for proteins in human tissues based on immunohistochemisty using tissue micro arrays.
getHpa1 <- getHpa(id, hpadata = "hpaNormalTissue")
getHpa11 <-getHpa1[which( getHpa1$Tissue%in% "lung"),]
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



View(datafull_join_getHpa10)

datafull_join_getHpa1a <- datafull_join_getHpa11 |>
  filter(Reliability == "Approved",
         Level == "High") |>
  arrange(desc("TissuesScore_Lung")) |>
  head()



View(datafull_join_getHpa1a)

datafull_join_getHpa1b <- datafull_join_getHpa11 |>
  filter(Reliability == "Enhanced",
         Level == "High") |>
  arrange(desc("TissuesScore_Lung")) |>
  head()
View(datafull_join_getHpa1b)

datafull_join_getHpa1c <-  datafull_join_getHpa11[order(datafull_join_getHpa11$TissuesScore_Lung, decreasing = TRUE), ]

View(datafull_join_getHpa1c)

datafull_join_getHpa1d <-  tissuesDB[order(tissuesDB$TissuesScore_Lung, decreasing = TRUE), ]
View(datafull_join_getHpa1d)

#head(getHpa(id, hpadata = "hpaNormalTissue"))
getHpa2 <-getHpa(id, hpadata = "hpaSubcellularLoc")
colnames(getHpa2)[2] <- "SYMBOL"
datafull_join_getHpa22 <- dplyr::full_join(datafull_join_getHpa11, 
                                           getHpa2, by='SYMBOL')


datafull_join_getHpa222 <- datafull_join_getHpa22 |>
  filter(Reliability.x == c("Approved","Enhanced")) |>
  arrange(desc("TissuesScore")) |>
  head()
View(datafull_join_getHpa222)
#head(getHpa(id, hpadata = "rnaGeneCellLine"))
getHpa3 <-getHpa(id, hpadata = "rnaGeneCellLine")
## ----getHpa2, eval=FALSE------------------------------------------------------
#  getHpa(id, type = "details")

