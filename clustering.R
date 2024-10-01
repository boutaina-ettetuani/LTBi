
###################################################################################################################
#####################################################################################################################
######################################################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GOSemSim)
library(DOSE)

gene <- unique(c(lit_entrezid$ENTREZID, anno_subset_.54992$ENTREZID,anno_subset_41055$ENTREZID, anno_subset_59184$ENTREZID))
gene_symbol_mat <- unique(c(lit_entrezid$SYMBOL, anno_subset_54992$SYMBOL,anno_subset_41055$SYMBOL, anno_subset_59184$SYMBOL))

ego <- enrichGO(gene  = gene,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                readable      = TRUE)

ego_trans <- enrichGO(gene  =  unique(c( anno_subset_54992$ENTREZID,anno_subset_41055$ENTREZID, anno_subset_59184$ENTREZID)),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                readable      = TRUE)

ego_lit <- enrichGO(gene  = lit_entrezid$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                readable      = TRUE)


DT::datatable(as.data.frame(ego@result[["Description"]]), options = list(scrollX = TRUE))
DT::datatable(as.data.frame(ego@result), options = list(scrollX = TRUE))

d <- godata('org.Hs.eg.db', ont=c("BP", "CC", "MF"))
ego2 <- pairwise_termsim(ego, method = "Wang", semData = d)
ego2_trans <- pairwise_termsim(ego_trans, method = "Wang", semData = d)
ego2_lit <- pairwise_termsim(ego_lit, method = "Wang", semData = d)
DT::datatable(as.data.frame(ego2@result[["Description"]]), options = list(scrollX = TRUE))

treeplot(ego2, showCategory = 30)
# use `hilight = FALSE` to remove ggtree::geom_hilight() layer.
treeplot(ego2, showCategory = 30, hilight = FALSE)
# use `offset` parameter to adjust the distance of bar and tree.
treeplot(ego2, showCategory = 30, hilight = FALSE, offset = 8)
# use `offset_tiplab` parameter to adjust the distance of nodes and branches.
treeplot(ego2, showCategory = 30, hilight = FALSE, offset_tiplab = 0.3)
keep <- rownames(ego2@termsim)[c(1:10, 16:20)]
keep
treeplot(ego2, showCategory = keep)
treeplot(ego2, showCategory = 20, 
         group_color = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442"))
########################################################
p2 <- treeplot(ego2, hclust_method = "average")

emapplot_cluster(ego2,showCategory = 200)




library(simplifyEnrichment)
set.seed(888)
go_id = random_GO(100)
#mat = GO_similarity(go_id, measure = "Wang")

mat = GO_similarity(ego2@result$ID, measure = "Wang")
mat_trans = GO_similarity(ego_trans@result$ID, measure = "Wang")
mat_lit = GO_similarity(ego2_lit@result$ID, measure = "Wang")
#df = simplifyGO(mat)
mat2 = GO_similarity(go_id, measure = "Wang")
library(simplifyEnrichment)
simplifyGO1 <- simplifyGO(mat)
simplifyGO_trans <- simplifyGO(mat_trans)
simplifyGO_lit <- simplifyGO(mat_lit)
simplifyGO(mat2)
cl = binary_cut(mat)
cl_trans = binary_cut(mat_trans)
cl_lit = binary_cut(mat_lit)

cl2 = binary_cut(mat2)
#export_to_shiny_app(mat, cl)

ht_clusters <- ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))
ht_clusters_trans <- ht_clusters(mat_trans, cl_trans, word_cloud_grob_param = list(max_width = 80))
ht_clusters_lit <- ht_clusters(mat_lit, cl_lit, word_cloud_grob_param = list(max_width = 80))

ht_clusters(mat2, cl2, word_cloud_grob_param = list(max_width = 80),
            order_by_size = TRUE)

compare_clustering_methods(mat2)
compare_clustering_methods(mat,method = "kmeans", plot_type = "heatmap")

df_mat = simplifyGO(mat, method = "kmeans")
df_trans = simplifyGO(mat_trans, method = "kmeans")
df_lit = simplifyGO(mat_lit, method = "kmeans")
############################################################################################
ego_mat <- as.data.frame(ego2@result)
df_mat_df <- as.data.frame(df_mat)
#############################################################################################
library(GOSemSim)
library(org.Hs.eg.db)
hsGO_all <- godata('org.Hs.eg.db', ont=c("BP", "CC", "MF"))
#############################################################################################
#############################################################################################
#morphogenesis,development
#############################################################################################
df_mat_1 <- subset(df_mat_df, df_mat_df$cluster == 1)
df_mat_1_subset <- ego_mat[which( ego_mat$ID%in% df_mat_1$id),]
dim(df_mat_1_subset)
head(df_mat_1_subset)
df_mat_1_subset_genes <-unlist(strsplit(df_mat_1_subset$geneID, split = "/"))
length(df_mat_1_subset_genes)
length(unique(df_mat_1_subset_genes))
head(df_mat_1_subset_genes)

##Expression
library(tidyselect)
df_mat_1_subset_genes_expression <- intersect(unique(df_mat_1_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_1_subset_genes_expression))
DataExpressioncsv_1_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_1_subset_genes_expression) )
dim(unique(DataExpressioncsv_1_subset_genes_cibled))
#############################################################################################
#############################################################################################
#development,morphogenesis
#############################################################################################
df_mat_2 <- subset(df_mat_df, df_mat_df$cluster == 2)
df_mat_2_subset <- ego_mat[which( ego_mat$ID%in% df_mat_2$id),]
dim(df_mat_2_subset)
head(df_mat_2_subset)
df_mat_2_subset_genes <-unlist(strsplit(df_mat_2_subset$geneID, split = "/"))
length(df_mat_2_subset_genes)
length(unique(df_mat_2_subset_genes))
head(df_mat_2_subset_genes)

##Expression

df_mat_2_subset_genes_expression <- intersect(unique(df_mat_2_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_2_subset_genes_expression))
DataExpressioncsv_2_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_2_subset_genes_expression) )
dim(unique(DataExpressioncsv_2_subset_genes_cibled))
#############################################################################################
#############################################################################################
#differentiation,cell,development
#############################################################################################
df_mat_3 <- subset(df_mat_df, df_mat_df$cluster == 3)
df_mat_3_subset <- ego_mat[which( ego_mat$ID%in% df_mat_3$id),]
dim(df_mat_3_subset)
head(df_mat_3_subset)
df_mat_3_subset_genes <-unlist(strsplit(df_mat_3_subset$geneID, split = "/"))
length(df_mat_3_subset_genes)
length(unique(df_mat_3_subset_genes))
head(df_mat_3_subset_genes)

##Expression

df_mat_3_subset_genes_expression <- intersect(unique(df_mat_3_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_3_subset_genes_expression))
DataExpressioncsv_3_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_3_subset_genes_expression) )
dim(unique(DataExpressioncsv_3_subset_genes_cibled))
#############################################################################################
#############################################################################################
#regulation,differentiation,cell,development
#############################################################################################
df_mat_4 <- subset(df_mat_df, df_mat_df$cluster == 4)
df_mat_4_subset <- ego_mat[which( ego_mat$ID%in% df_mat_4$id),]
dim(df_mat_4_subset)
head(df_mat_4_subset)
df_mat_4_subset_genes <-unlist(strsplit(df_mat_4_subset$geneID, split = "/"))
length(df_mat_4_subset_genes)
length(unique(df_mat_4_subset_genes))
head(df_mat_4_subset_genes)

##Expression

df_mat_4_subset_genes_expression <- intersect(unique(df_mat_4_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_4_subset_genes_expression))
DataExpressioncsv_4_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_4_subset_genes_expression) )
dim(unique(DataExpressioncsv_4_subset_genes_cibled))

#############################################################################################
#############################################################################################
#production, regulation
#############################################################################################
df_mat_5 <- subset(df_mat_df, df_mat_df$cluster == 5)
df_mat_5_subset <- ego_mat[which( ego_mat$ID%in% df_mat_5$id),]
dim(df_mat_5_subset)
head(df_mat_5_subset)
df_mat_5_subset_genes <-unlist(strsplit(df_mat_5_subset$geneID, split = "/"))
length(df_mat_5_subset_genes)
length(unique(df_mat_5_subset_genes))
head(df_mat_5_subset_genes)

##Expression

df_mat_5_subset_genes_expression <- intersect(unique(df_mat_5_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_5_subset_genes_expression))
DataExpressioncsv_5_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_5_subset_genes_expression) )
dim(unique(DataExpressioncsv_5_subset_genes_cibled))
#############################################################################################
#############################################################################################
#regulation,proteine,activity,process,positive,proteine, negative
#############################################################################################
df_mat_6 <- subset(df_mat_df, df_mat_df$cluster == 6)
df_mat_6_subset <- ego_mat[which( ego_mat$ID%in% df_mat_6$id),]
dim(df_mat_6_subset)
head(df_mat_6_subset)
df_mat_6_subset_genes <-unlist(strsplit(df_mat_6_subset$geneID, split = "/"))
length(df_mat_6_subset_genes)
length(unique(df_mat_6_subset_genes))
head(df_mat_6_subset_genes)

##Expression

df_mat_6_subset_genes_expression <- intersect(unique(df_mat_6_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_6_subset_genes_expression))
DataExpressioncsv_6_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_6_subset_genes_expression) )
dim(unique(DataExpressioncsv_6_subset_genes_cibled))
#############################################################################################
#############################################################################################
#s
#############################################################################################
df_mat_7 <- subset(df_mat_df, df_mat_df$cluster == 7)
df_mat_7_subset <- ego_mat[which( ego_mat$ID%in% df_mat_7$id),]
dim(df_mat_7_subset)
head(df_mat_7_subset)
df_mat_7_subset_genes <-unlist(strsplit(df_mat_7_subset$geneID, split = "/"))
length(df_mat_7_subset_genes)
length(unique(df_mat_7_subset_genes))
head(df_mat_7_subset_genes)

##Expression

df_mat_7_subset_genes_expression <- intersect(unique(df_mat_7_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_7_subset_genes_expression))
DataExpressioncsv_7_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_7_subset_genes_expression) )
dim(unique(DataExpressioncsv_7_subset_genes_cibled))
#############################################################################################
#############################################################################################
#protein,modification,histone,metabolic,peptidyl,methylation
#############################################################################################
df_mat_8 <- subset(df_mat_df, df_mat_df$cluster == 8)
df_mat_8_subset <- ego_mat[which( ego_mat$ID%in% df_mat_8$id),]
dim(df_mat_8_subset)
head(df_mat_8_subset)
df_mat_8_subset_genes <-unlist(strsplit(df_mat_8_subset$geneID, split = "/"))
length(df_mat_8_subset_genes)
length(unique(df_mat_8_subset_genes))
head(df_mat_8_subset_genes)

##Expression

df_mat_8_subset_genes_expression <- intersect(unique(df_mat_8_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_8_subset_genes_expression))
DataExpressioncsv_8_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_8_subset_genes_expression) )
dim(unique(DataExpressioncsv_8_subset_genes_cibled))
#############################################################################################
#############################################################################################
#process,metabolic,biosynthetic
#############################################################################################
df_mat_9 <- subset(df_mat_df, df_mat_df$cluster == 9)
df_mat_9_subset <- ego_mat[which( ego_mat$ID%in% df_mat_9$id),]
dim(df_mat_9_subset)
head(df_mat_9_subset)
df_mat_9_subset_genes <-unlist(strsplit(df_mat_9_subset$geneID, split = "/"))
length(df_mat_9_subset_genes)
length(unique(df_mat_9_subset_genes))
head(df_mat_9_subset_genes)

##Expression

df_mat_9_subset_genes_expression <- intersect(unique(df_mat_9_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_9_subset_genes_expression))
DataExpressioncsv_9_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_9_subset_genes_expression) )
dim(unique(DataExpressioncsv_9_subset_genes_cibled))
#############################################################################################
#############################################################################################
#process,biosynthetic,metabolic,catabolic
#############################################################################################
df_mat_10 <- subset(df_mat_df, df_mat_df$cluster == 10)
df_mat_10_subset <- ego_mat[which( ego_mat$ID%in% df_mat_10$id),]
dim(df_mat_10_subset)
head(df_mat_10_subset)
df_mat_10_subset_genes <-unlist(strsplit(df_mat_10_subset$geneID, split = "/"))
length(df_mat_10_subset_genes)
length(unique(df_mat_10_subset_genes))
head(df_mat_10_subset_genes)

##Expression

df_mat_10_subset_genes_expression <- intersect(unique(df_mat_10_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_10_subset_genes_expression))
DataExpressioncsv_10_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_10_subset_genes_expression) )
dim(unique(DataExpressioncsv_10_subset_genes_cibled))
#############################################################################################
#############################################################################################
#process,metabolic,biosynthetic
#############################################################################################
df_mat_11 <- subset(df_mat_df, df_mat_df$cluster == 11)
df_mat_11_subset <- ego_mat[which( ego_mat$ID%in% df_mat_11$id),]
dim(df_mat_11_subset)
head(df_mat_11_subset)
df_mat_11_subset_genes <-unlist(strsplit(df_mat_11_subset$geneID, split = "/"))
length(df_mat_11_subset_genes)
length(unique(df_mat_11_subset_genes))
head(df_mat_11_subset_genes)

##Expression

df_mat_11_subset_genes_expression <- intersect(unique(df_mat_11_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_11_subset_genes_expression))
DataExpressioncsv_11_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_11_subset_genes_expression) )
dim(unique(DataExpressioncsv_11_subset_genes_cibled))
#############################################################################################
#############################################################################################
#signaling,regulation,pathway,receptor,reponse,protein
#############################################################################################
df_mat_12 <- subset(df_mat_df, df_mat_df$cluster == 12)
df_mat_12_subset <- ego_mat[which( ego_mat$ID%in% df_mat_12$id),]
dim(df_mat_12_subset)
head(df_mat_12_subset)
df_mat_12_subset_genes <-unlist(strsplit(df_mat_12_subset$geneID, split = "/"))
length(df_mat_12_subset_genes)
length(unique(df_mat_12_subset_genes))
head(df_mat_12_subset_genes)

##Expression

df_mat_12_subset_genes_expression <- intersect(unique(df_mat_12_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_12_subset_genes_expression))
DataExpressioncsv_12_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_12_subset_genes_expression) )
dim(unique(DataExpressioncsv_12_subset_genes_cibled))
#############################################################################################
#############################################################################################
#regulation,transport,protein,transmembrane,membrane,secrection,localization
#############################################################################################
df_mat_13 <- subset(df_mat_df, df_mat_df$cluster == 13)
df_mat_13_subset <- ego_mat[which( ego_mat$ID%in% df_mat_13$id),]
dim(df_mat_13_subset)
head(df_mat_13_subset)
df_mat_13_subset_genes <-unlist(strsplit(df_mat_13_subset$geneID, split = "/"))
length(df_mat_13_subset_genes)
length(unique(df_mat_13_subset_genes))
head(df_mat_13_subset_genes)

##Expression

df_mat_13_subset_genes_expression <- intersect(unique(df_mat_13_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_13_subset_genes_expression))
DataExpressioncsv_13_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_13_subset_genes_expression) )
dim(unique(DataExpressioncsv_13_subset_genes_cibled))
#############################################################################################
#############################################################################################
#reponse,cellular,stimulus
#############################################################################################
df_mat_14 <- subset(df_mat_df, df_mat_df$cluster == 14)
df_mat_14_subset <- ego_mat[which( ego_mat$ID%in% df_mat_14$id),]
dim(df_mat_14_subset)
head(df_mat_14_subset)
df_mat_14_subset_genes <-unlist(strsplit(df_mat_14_subset$geneID, split = "/"))
length(df_mat_14_subset_genes)
length(unique(df_mat_14_subset_genes))
head(df_mat_14_subset_genes)

##Expression

df_mat_14_subset_genes_expression <- intersect(unique(df_mat_14_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_14_subset_genes_expression))
DataExpressioncsv_14_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_14_subset_genes_expression) )
dim(unique(DataExpressioncsv_14_subset_genes_cibled))
#############################################################################################
#############################################################################################
#regulation cell,immune,activation,response
#############################################################################################
df_mat_15 <- subset(df_mat_df, df_mat_df$cluster == 15)
df_mat_15_subset <- ego_mat[which( ego_mat$ID%in% df_mat_15$id),]
dim(df_mat_15_subset)
head(df_mat_15_subset)
df_mat_15_subset_genes <-unlist(strsplit(df_mat_15_subset$geneID, split = "/"))
length(df_mat_15_subset_genes)
length(unique(df_mat_15_subset_genes))
head(df_mat_15_subset_genes)

##Expression

df_mat_15_subset_genes_expression <- intersect(unique(df_mat_15_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_15_subset_genes_expression))
DataExpressioncsv_15_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_15_subset_genes_expression) )
dim(unique(DataExpressioncsv_15_subset_genes_cibled))
#############################################################################################
#############################################################################################
#regulation,muscle,contraction,cardiac,blood,cardiac
#############################################################################################
df_mat_16 <- subset(df_mat_df, df_mat_df$cluster == 16)
df_mat_16_subset <- ego_mat[which( ego_mat$ID%in% df_mat_16$id),]
dim(df_mat_16_subset)
head(df_mat_16_subset)
df_mat_16_subset_genes <-unlist(strsplit(df_mat_16_subset$geneID, split = "/"))
length(df_mat_16_subset_genes)
length(unique(df_mat_16_subset_genes))
head(df_mat_16_subset_genes)

##Expression

df_mat_16_subset_genes_expression <- intersect(unique(df_mat_16_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_16_subset_genes_expression))
DataExpressioncsv_16_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_16_subset_genes_expression) )
dim(unique(DataExpressioncsv_16_subset_genes_cibled))
#############################################################################################
#############################################################################################
#organization,membrane,microtubule,mitochondrial,regulation,protein,cell
#############################################################################################
df_mat_17 <- subset(df_mat_df, df_mat_df$cluster == 17)
df_mat_17_subset <- ego_mat[which( ego_mat$ID%in% df_mat_17$id),]
dim(df_mat_17_subset)
head(df_mat_17_subset)
df_mat_17_subset_genes <-unlist(strsplit(df_mat_17_subset$geneID, split = "/"))
length(df_mat_17_subset_genes)
length(unique(df_mat_17_subset_genes))
head(df_mat_17_subset_genes)

##Expression

df_mat_17_subset_genes_expression <- intersect(unique(df_mat_17_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_17_subset_genes_expression))
DataExpressioncsv_17_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_17_subset_genes_expression) )
dim(unique(DataExpressioncsv_17_subset_genes_cibled))
#############################################################################################
#############################################################################################
#cell,process,localization,poliferation,adhesion,cellular,migration
#############################################################################################
df_mat_18 <- subset(df_mat_df, df_mat_df$cluster == 18)
df_mat_18_subset <- ego_mat[which( ego_mat$ID%in% df_mat_18$id),]
dim(df_mat_18_subset)
head(df_mat_18_subset)
df_mat_18_subset_genes <-unlist(strsplit(df_mat_18_subset$geneID, split = "/"))
length(df_mat_18_subset_genes)
length(unique(df_mat_18_subset_genes))
head(df_mat_18_subset_genes)

##Expression

df_mat_18_subset_genes_expression <- intersect(unique(df_mat_18_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_18_subset_genes_expression))
DataExpressioncsv_18_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_18_subset_genes_expression) )
dim(unique(DataExpressioncsv_18_subset_genes_cibled))
#############################################################################################
#############################################################################################
#regulation,organization,assembly
#############################################################################################
df_mat_19 <- subset(df_mat_df, df_mat_df$cluster == 19)
df_mat_19_subset <- ego_mat[which( ego_mat$ID%in% df_mat_19$id),]
dim(df_mat_19_subset)
head(df_mat_19_subset)
df_mat_19_subset_genes <-unlist(strsplit(df_mat_19_subset$geneID, split = "/"))
length(df_mat_19_subset_genes)
length(unique(df_mat_19_subset_genes))
head(df_mat_19_subset_genes)

##Expression

df_mat_19_subset_genes_expression <- intersect(unique(df_mat_19_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_19_subset_genes_expression))
DataExpressioncsv_19_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_19_subset_genes_expression) )
dim(unique(DataExpressioncsv_19_subset_genes_cibled))


#############################################################################################
#############################################################################################
#regulation,cell,positive
#############################################################################################
df_mat_20 <- subset(df_mat_df, df_mat_df$cluster == 20)
df_mat_20_subset <- ego_mat[which( ego_mat$ID%in% df_mat_20$id),]
dim(df_mat_20_subset)
head(df_mat_20_subset)
df_mat_20_subset_genes <-unlist(strsplit(df_mat_20_subset$geneID, split = "/"))
length(df_mat_20_subset_genes)
length(unique(df_mat_20_subset_genes))
head(df_mat_20_subset_genes)

##Expression

df_mat_20_subset_genes_expression <- intersect(unique(df_mat_20_subset_genes),unique(gene_symbol_mat))
length(unique(df_mat_20_subset_genes_expression))
DataExpressioncsv_20_subset_genes_cibled <- subset(datafull_join,datafull_join$SYMBOL %in% unique(df_mat_20_subset_genes_expression) )
dim(unique(DataExpressioncsv_20_subset_genes_cibled))




tst <- c(unique(df_mat_1_subset_genes),unique(gene_symbol_mat))
tst <- tst[duplicated(tst)]
tst[duplicated(tst)]







library(VennDiagram)
library(purrr)
library(RVenn)
library(ggplot2)
toy = list(dt4= dt4_subset_annot_genes$SYMBOL,dt5 = dt5_subset_annot_genes$SYMBOL, dt10 =dt10_subset_annot_genes$SYMBOL, dt9 =dt9_subset_annot_genes$SYMBOL)
toy = Venn(toy)

#####Draw the Venn Diagram
ggvenn(toy, slice = c(1, 2,3))
ggvenn(toy, slice = c(1, 2))
ggvenn(toy, slice = c(1,3))

ggvenn(toy, slice = c(1,4))
ggvenn(toy, slice = c(2,4))
ggvenn(toy, slice = c(3,4))
ggvenn(toy, slice = c(2,1))
ggvenn(toy, slice = c(1, 2, 3),fill = c("deeppink", "dodgerblue3", "pink"))
#fill = c("gold", "dodgerblue3", "deeppink")
###################################################################################################################
###################################################################################################################
datafull_join <- dplyr::full_join(dt4_subset_annot_genes[,c("SYMBOL","score")], 
                                  dt5_subset_annot_genes[,c("SYMBOL","score")], by='SYMBOL') %>%
  dplyr::full_join(., dt9_subset_annot_genes[,c("SYMBOL","score")], by='SYMBOL') %>%
  dplyr::full_join(., dt10_subset_annot_genes[,c("SYMBOL","score")], by='SYMBOL') 

colnames(datafull_join)[2] <- "direct_CL_4"
colnames(datafull_join)[3] <- "direct_CL_5"
colnames(datafull_join)[4] <- "Undirect_CL_9"
colnames(datafull_join)[5] <- "Undirect_CL_10"

head(datafull_join)
library("plyr") 

dfr <- rbind(dt4_subset_annot_genes[,c("SYMBOL","score","expression")],dt5_subset_annot_genes[,c("SYMBOL","score","expression")],
             dt9_subset_annot_genes[,c("SYMBOL","score","expression")],dt10_subset_annot_genes[,c("SYMBOL","score","expression")])


colnames(dfr)[3] <- "Cluster"
#######################3
require(FactoMineR)
require(factoextra)
require(ggplot2)
require(tidyr)
require(dplyr)
require(MASS)
require(reshape2)
require(cowplot)
# Run the PCA
datadfr_dt4_subset_annot_genes=dt4_subset_annot_genes[,c("score","expression")]
row.names(datadfr_dt4_subset_annot_genes)=dt4_subset_annot_genes$SYMBOL

datadfr_dt5_subset_annot_genes=dt5_subset_annot_genes[,c("score","expression")]
row.names(datadfr_dt5_subset_annot_genes)=dt5_subset_annot_genes$SYMBOL

datadfr_dt9_subset_annot_genes=dt9_subset_annot_genes[,c("score","expression")]
row.names(datadfr_dt9_subset_annot_genes)=dt9_subset_annot_genes$SYMBOL

datadfr_dt10_subset_annot_genes=dt10_subset_annot_genes[,c("score","expression")]
row.names(datadfr_dt10_subset_annot_genes)=dt10_subset_annot_genes$SYMBOL



datadfr_dt4_subset_annot_genes$clustId <- "direct_CL_4"
datadfr_dt5_subset_annot_genes$clustId <- "direct_CL_5"
datadfr_dt9_subset_annot_genes$clustId <- "Undirect_CL_9"
datadfr_dt10_subset_annot_genes$clustId <- "Undirect_CL_10"


rbind1 <-  rbind(datadfr_dt4_subset_annot_genes[,c("score","clustId")],datadfr_dt5_subset_annot_genes[,c("score","clustId")],datadfr_dt9_subset_annot_genes[,c("score","clustId")],datadfr_dt10_subset_annot_genes[,c("score","clustId")])
