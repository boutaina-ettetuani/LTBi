
library(GOxploreR)
# Retrieve the level of a GO biological process term
#regulation
goterms0 <- unique(c(df_mat_4_subset$ID,df_mat_6_subset$ID,df_mat_7_subset$ID,df_mat_13_subset$ID,df_mat_15_subset$ID,df_mat_16_subset$ID,
                    df_mat_17_subset$ID,df_mat_18_subset$ID))
df_mat_Regu_subset0 <- rbind(df_mat_4_subset,df_mat_6_subset,df_mat_7_subset,df_mat_13_subset,df_mat_15_subset,
                            df_mat_16_subset,df_mat_17_subset,df_mat_19_subset,df_mat_20_subset)

goterms <- unique(c(df_mat_4_subset$ID,df_mat_6_subset$ID,df_mat_7_subset$ID,df_mat_13_subset$ID,df_mat_15_subset$ID,
                     df_mat_17_subset$ID,df_mat_19_subset$ID))
#################################################################################################################
df_mat_Regu_subset0 <- rbind(df_mat_4_subset,df_mat_6_subset,df_mat_7_subset,df_mat_13_subset,df_mat_15_subset,
                             df_mat_17_subset,df_mat_19_subset)
df_mat_Regu_subset <- rbind(df_mat_4_subset,df_mat_5_subset,df_mat_6_subset,df_mat_7_subset,df_mat_12_subset,df_mat_13_subset,df_mat_15_subset,
                            df_mat_16_subset,df_mat_17_subset,df_mat_19_subset,df_mat_20_subset)
#####################################################
length(unique(df_mat_Regu_subset$ID))
df_mat_Regu_gene_subset <- unlist(strsplit(df_mat_Regu_subset$geneID, split = "/"))
length(df_mat_Regu_gene_subset)
length(unique(df_mat_Regu_gene_subset))
length(unique(goterms))

df_mat_Regu_subset_genes <-unlist(strsplit(df_mat_Regu_subset$geneID, split = "/"))
length(df_mat_Regu_subset_genes)
length(unique(df_mat_Regu_subset_genes))
head(df_mat_Regu_subset_genes)

#other
#goterms <- unique(unique(df_mat_Regu_subset$ID))
#length(goterms)
goterms <- unique(c(df_mat_4_subset$ID,df_mat_6_subset$ID,df_mat_7_subset$ID,df_mat_13_subset$ID,df_mat_15_subset$ID,
                    df_mat_17_subset$ID,df_mat_19_subset$ID))

length(goterms)
#write.csv(goterms, file = 'C:/Users/hp/Desktop/ltbiTranscrp/goterms.csv', row.names = FALSE)

GOTermBPOnLevel_reg <- GOTermBPOnLevel(goterm = unique(goterms))
#set.seed (9996743)
#4
GOTermBPOnLevel_reg_4 <- GOTermBPOnLevel(goterm = unique(df_mat_4_subset$ID))
GOTermBPOnLevel_reg_4
length(unique(GOTermBPOnLevel_reg$Term))
# Returns the categories of the GO-terms in the list
getGOcategory_reg <- getGOcategory(goterm = goterms)
getGOcategory_reg
visRsubDAGBP_reg <- visRsubDAGBP(goterm = goterms, organism = "Human")
visRsubDAGBP_reg

#visRsubDAGBP(goterm = df3_4_subset_cibled$ID, organism = "Human")
#visRsubDAGBP(goterm = df3_10_subset_cibled$ID[1:24], organism = "Human")


# Ordering of the GO-terms in the list
distRankingGO(goterm = goterms, domain = "BP", plot = TRUE)



GOTermBPOnLevel
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 4)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 5)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 6)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 7)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 8)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 9)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 10)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 11)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 12)
filter(GOTermBPOnLevel, GOTermBPOnLevel$Level == 13)

# Ordering of the GO-terms in the list
scoreRankingGO <- scoreRankingGO(goterm = goterms, domain = "BP", plot = FALSE)
scoreRankingGO
# We Prioritize the given biological process GO-terms
prioritizedGOTerms <- prioritizedGOTerms(lst = goterms, organism = "Human", sp = TRUE, domain = "BP")


library(systemPipeR)


#GOTermBPOnLevel =unique(GOTermBPOnLevel$Term
#df_mat_8_subset
length(unique(df_mat_Regu_gene_subset))
length(unique(goterms))
length(unique(df_mat_Regu_subset$ID))
length(unique(GOTermBPOnLevel$Term))
setlist_symbol <- list(goterms=goterms,df_mat_Regu_subset =unique(df_mat_Regu_subset$ID))
setlist_symbol_go <- list(GOTermBPOnLevel =unique(GOTermBPOnLevel_reg$Term),df_mat_Regu_subset =unique(df_mat_Regu_subset$ID))

## 2-way Venn diagram
vennset2_symbol_go <- overLapper(setlist_symbol_go[1:2], type="vennsets")
#vennset2_symbol <- overLapper(setlist_symbol[c(1,3)], type="vennsets")

vennPlot(vennset2_symbol_go)

vennset2_symbol_go_venn = Venn(setlist_symbol_go[1:2])
overlap(vennset2_symbol_go_venn)
length(unique(overlap(vennset2_symbol_go_venn)))



df_mat_Regu_subset_overlap <- df_mat_Regu_subset[which( df_mat_Regu_subset$ID%in% unique(overlap(vennset2_symbol_go_venn))),]
dim(df_mat_Regu_subset_overlap)
head(df_mat_Regu_subset_overlap)
df_mat_Regu_subset_overlap_genes <-unlist(strsplit(df_mat_Regu_subset_overlap$geneID, split = "/"))
length(unique(df_mat_Regu_subset_overlap_genes))




toy_ensembl <-  AnnotationDbi::select(hgu133plus2.db,
                                      keys = unique(df_mat_Regu_subset_genes),#unique(df_mat_Regu_subset_overlap_genes)
                                      columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                      keytype = "SYMBOL")

library(hpar)
id <- unique(toy_ensembl$ENSEMBL)
getHpa <- getHpa(id, hpadata = "hpaNormalTissue")

getHpa

library(disgenet2r)
library(RCurl)
library(httr)
disgenet_api_key <- get_disgenet_api_key(
  email = "ettetuani.boutaina@gmail.com", 
  password = "Boutaina@ett1992" )
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

length(unique(toy_ensembl$ENTREZID))
data1 <- gene2disease( gene = unique(toy_ensembl$ENTREZID[1:900]), vocabulary = "ENTREZ",
                       database = "CURATED")
data1 <- gene2disease( gene = unique(toy_ensembl$SYMBOL[1:750]), verbose = TRUE)

results <- extract(data1)
##################################################################################################
data2 <- gene2disease( gene = toy_ensembl$SYMBOL[1:750], verbose = TRUE)
results2 <- extract(data2)
plot( data2,
      class = "Network",
      prop = 20)
plot( data2,
      class = "DiseaseClass",
      prop = 3)

data3 <- gene2disease( gene = toy_ensembl$SYMBOL[750:1472], verbose = TRUE)
results3 <- extract(data3)

data4 <- gene2disease( gene = toy_ensembl$SYMBOL[1472:2000], verbose = TRUE)
results4 <- extract(data4)

data5 <- gene2disease( gene = toy_ensembl$SYMBOL[2000:2700], verbose = TRUE)
results5 <- extract(data5)

data6 <- gene2disease( gene = toy_ensembl$SYMBOL[2700:3400], verbose = TRUE)
results6 <- extract(data6)

data7 <- gene2disease( gene = toy_ensembl$SYMBOL[3400:4100], verbose = TRUE)
results7 <- extract(data7)

data8 <- gene2disease( gene = toy_ensembl$SYMBOL[4100:4800], verbose = TRUE)
results8 <- extract(data8)

data9 <- gene2disease( gene = toy_ensembl$SYMBOL[4800:5500], verbose = TRUE)
results9 <- extract(data9)

data10 <- gene2disease( gene = toy_ensembl$SYMBOL[5500:6200], verbose = TRUE)
results10 <- extract(data10)

data11 <- gene2disease( gene = toy_ensembl$SYMBOL[6200:6800], verbose = TRUE)
results11 <- extract(data11)

data12 <- gene2disease( gene = toy_ensembl$SYMBOL[6800:7200], verbose = TRUE)
results12 <- extract(data12)

data13 <- gene2disease( gene = toy_ensembl$SYMBOL[7200:7610], verbose = TRUE)
results13 <- extract(data13)

data_reg_gene2disease_data = cbind(data2,data3,data4,data5,data6,data7,data8,data9,data10,data11,data12,data13)

data_reg_gene2disease = rbind(results2,results3,results4,results5,results6,results7,results8,results9,results10,results11,results12,results13)
View(data_reg_gene2disease)


data_reg_gene2disease_TB_TBi <- rbind(data_reg_gene2disease_TB,data_reg_gene2disease_TBi)

data_reg_gene2disease_TB <- dplyr::filter(data_reg_gene2disease, grepl("Tuberculosis", disease_name, ignore.case = TRUE))
data_reg_gene2disease_TBi <- dplyr::filter(data_reg_gene2disease_TB, grepl("infections", disease_class_name, ignore.case = TRUE))
###plot
plot( data_reg_gene2disease_TBi,
      class = "DiseaseClass",
      prop = 3)
# View the new data set
head(data_reg_gene2disease_TBi)
#https://www.disgenet.org/static/disgenet2r/disgenet2r.html
#write.csv(data_reg_gene2disease_TBi, file = 'C:/Users/hp/Desktop/ltbiTranscrp/data_reg_gene2disease_TBi.csv', row.names = FALSE)

#YES
#7610#####################################################################################################

data22 <- gene2evidence( gene = toy_ensembl$SYMBOL[1:20],
                        vocabulary = "HGNC",
                        database = "ALL")
results22 <- extract(data22)
#NO
#####################################################################################################
#####################################################################################################
#####################################################################################################



x <- enrichDO(gene          = unique(toy_ensembl$ENTREZID),
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
#View(x@result)

#NO
#######################################################################################################3
#YES
dgn <- enrichDGN(unique(toy_ensembl$ENTREZID)) 
head(dgn)
View(dgn@result)


enrichDGN_TB_pulmo <- dplyr::filter(dgn@result, grepl("Tuberculosis, Pulmonary", Description, ignore.case = TRUE))
# View the new data set
head(enrichDGN_TB_pulmo)
#write.csv(enrichDGN_TB_pulmo, file = 'C:/Users/hp/Desktop/ltbiTranscrp/enrichDGN_TB_pulmo.csv', row.names = FALSE)

enrichDGN_LTBi <- dplyr::filter(dgn@result, grepl("Latent Tuberculosis", Description, ignore.case = TRUE))
# View the new data set
head(enrichDGN_LTBi)

enrichDGN_TB_active <- dplyr::filter(dgn@result, grepl("Active Tuberculosis", Description, ignore.case = TRUE))
# View the new data set
head(enrichDGN_TB_active)
#write.csv(enrichDGN_LTBi, file = 'C:/Users/hp/Desktop/ltbiTranscrp/enrichDGN_LTBi.csv', row.names = FALSE)

#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
##############plot 
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
enrichDGN_LTBi_plot <- dplyr::filter(dgn, grepl("Tuberculosis", dgn@result$Description, ignore.case = TRUE))
##plot1
library(enrichplot)
barplot(enrichDGN_LTBi_plot, showCategory=40)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr) 
#mutate(enrichDGN_LTBi_plot, qscore = -log(p.adjust, base=10)) %>% barplot(x="qscore")
## convert gene ID to Symbol
##plot2
edox <- setReadable(enrichDGN_LTBi_plot, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue")
p3 <- cnetplot(edox, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
##plot3
p11 <- cnetplot(edox, node_label="category", 
               cex_label_category = 1.2) 
p22 <- cnetplot(edox, node_label="gene", 
               cex_label_gene = 0.8) 
p33 <- cnetplot(edox, node_label="all") 
p44 <- cnetplot(edox, node_label="none", 
               color_category='firebrick', 
               color_gene='steelblue') 
cowplot::plot_grid(p11, p22, p33, p44, ncol=2, labels=LETTERS[1:4])
##plot4
p111 <- heatplot(edox, showCategory=5)
p222 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p111, p222, ncol=1, labels=LETTERS[1:2])
###plot5
edox2 <- pairwise_termsim(edox)
p1_1 <- treeplot(edox2)
p2_2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1_1, p2_2, tag_levels='A')
################
###end of the plot
# View some of the values
#################################################enrichDGN_TB_pulmo
head(enrichDGN_TB_pulmo)
Entrez_genes_enrichDGN_TB_pulmo <- unlist(strsplit(enrichDGN_TB_pulmo$geneID, split = "/"))
unique(Entrez_genes_enrichDGN_TB_pulmo)
GENENAME_genes_enrichDGN_TB_pulmo <-  AnnotationDbi::select(hgu133plus2.db,
                                                        keys = unique(Entrez_genes_enrichDGN_TB_pulmo),#unique(df_mat_Regu_subset_overlap_genes)
                                                        columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                                        keytype = "ENTREZID")

unique(GENENAME_genes_enrichDGN_TB_pulmo)
unique(GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL)
#write.csv(GENENAME_genes_enrichDGN_TB_pulmo, file = 'C:/Users/hp/Desktop/ltbiTranscrp/GENENAME_genes_enrichDGN_TB_pulmo.csv', row.names = FALSE)


####################################################ltbi
head(enrichDGN_LTBi)
Entrez_genes_enrichDGN_LTBi <- unlist(strsplit(enrichDGN_LTBi$geneID, split = "/"))
unique(Entrez_genes_enrichDGN_LTBi)
GENENAME_genes_enrichDGN_LTBi <-  AnnotationDbi::select(hgu133plus2.db,
                                                        keys = unique(Entrez_genes_enrichDGN_LTBi),#unique(df_mat_Regu_subset_overlap_genes)
                                                        columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                                        keytype = "ENTREZID")

unique(GENENAME_genes_enrichDGN_LTBi)
unique(GENENAME_genes_enrichDGN_LTBi$SYMBOL)
#write.csv(GENENAME_genes_enrichDGN_LTBi, file = 'C:/Users/hp/Desktop/ltbiTranscrp/GENENAME_genes_enrichDGN_LTBi.csv', row.names = FALSE)
####################################################ltbi
head(enrichDGN_TB_active)
Entrez_genes_enrichDGN_TB_active <- unlist(strsplit(enrichDGN_TB_active$geneID, split = "/"))
unique(Entrez_genes_enrichDGN_TB_active)
GENENAME_genes_enrichDGN_TB_active <-  AnnotationDbi::select(hgu133plus2.db,
                                                        keys = unique(Entrez_genes_enrichDGN_TB_active),#unique(df_mat_Regu_subset_overlap_genes)
                                                        columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                                        keytype = "ENTREZID")

unique(GENENAME_genes_enrichDGN_TB_active)
unique(GENENAME_genes_enrichDGN_TB_active$SYMBOL)

####################################################################33
library(hpar)
id <- unique(GENENAME_genes_enrichDGN_LTBi$ENSEMBL)
getHpa <- getHpa(id, hpadata = "hpaNormalTissue")

getHpa
#NO
##########################################################################################################
library(clusterProfiler)
ns.kegg2 <- enrichKEGG(gene = unique(toy_ensembl$ENTREZID[1:750]),
                      organism = 'hsa',
                      pvalueCutoff = 0.05)

##############################################################################################################



library(VennDiagram)
library(purrr)
library(RVenn)
library(ggplot2)
toy = list(genetodisease_TBi= data_reg_gene2disease_TBi$gene_symbol,enrichment_LTBi = GENENAME_genes_enrichDGN_LTBi$SYMBOL, enrichment_Pulmonary_TB = GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL,enrichDGN_TB_active = GENENAME_genes_enrichDGN_TB_active$SYMBOL )
toy = Venn(toy)
overlap(toy)

toy2 <- overLapper(toy[1:3], type="vennsets")
vennPlot(toy2)
toy2@vennlist[["genetodisease_TBi_enrichment_LTBi_enrichment_Pulmonary_TB"]]
#####Draw the Venn Diagram
ggvenn(toy, slice = c(1, 2,3))

####3
library(systemPipeR)
list_LTBi <- list(genetodisease_TBi= unique(data_reg_gene2disease_TBi$gene_symbol),enrichment_LTBi = unique(GENENAME_genes_enrichDGN_LTBi$SYMBOL), enrichment_Pulmonary_TB = unique(GENENAME_genes_enrichDGN_TB_pulmo$SYMBOL),enrichDGN_TB_active = unique(GENENAME_genes_enrichDGN_TB_active$SYMBOL))
list_LTBi_entrez <- list(genetodisease_TBi= unique(data_reg_gene2disease_TBi$geneid),enrichment_LTBi = unique(GENENAME_genes_enrichDGN_LTBi$ENTREZID), enrichment_Pulmonary_TB = unique(GENENAME_genes_enrichDGN_TB_pulmo$ENTREZID),enrichDGN_TB_active = unique(GENENAME_genes_enrichDGN_TB_active$ENTREZID))

## 2-way Venn diagram
Results <- overLapper(list_LTBi[1:4], type="vennsets")
vennPlot(Results)

GENENAME_genes_enrichDGN_LTBi_43 <-GENENAME_genes_enrichDGN_LTBi[which( GENENAME_genes_enrichDGN_LTBi$SYMBOL%in% Results@vennlist[["enrichment_LTBi"]]),]
unique(GENENAME_genes_enrichDGN_LTBi_43)
unique(GENENAME_genes_enrichDGN_LTBi_43$SYMBOL)
write.csv(GENENAME_genes_enrichDGN_LTBi_43, file = 'C:/Users/hp/Desktop/ltbiTranscrp/GENENAME_genes_enrichDGN_LTBi_43.csv', row.names = FALSE)


gene2dis_data1 <- gene2disease( gene = c(list_LTBi$genetodisease_TBi,list_LTBi$enrichment_LTBi), verbose = TRUE)
results_gene2dis_data1 <- extract(gene2dis_data1)
gene2dis_data1@qresult
gene2dis_data2 <- gene2disease( gene = c(list_LTBi$enrichment_Pulmonary_TB,list_LTBi$enrichDGN_TB_active), verbose = TRUE)
results_gene2dis_data2 <- extract(gene2dis_data2)
results_gene2dis_data  <- rbind(gene2dis_data1@qresult,gene2dis_data2@qresult)
View(results_gene2dis_data)
length(unique(results_gene2dis_data$gene_symbol))

get_Genes_info= getBM(attributes='hgnc_symbol',
                   filters = c('hgnc_symbol','chromosome_name','chromosomal_region','start','end','strand'), 
                   values = unique(results_gene2dis_data$gene_symbol),
                   mart = mart,
                   useCache = FALSE) 

get_Genes_info= getBM(attributes='hgnc_symbol',
                      filters = c('hgnc_symbol','start','end'), 
                      values = unique(results_gene2dis_data$gene_symbol),
                      mart = mart,
                      useCache = FALSE) 
get_Genes_info

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
#listAttributes(ensembl)
chr11_genes <- getBM(attributes=c('hgnc_symbol','5utr'), filters =
                        'hgnc_symbol', values ="PPARG", mart = ensembl)
View(chr11_genes)
chr1_genes <- getBM(attributes=c('hgnc_symbol','chromosome_name','gene_exon_intron','band'), filters =
                      'hgnc_symbol', values =unique(results_gene2dis_data$gene_symbol[1:90]), mart = ensembl)
chr2_genes <- getBM(attributes=c('hgnc_symbol','chromosome_name','gene_exon_intron','band'), filters =
                      'hgnc_symbol', values =unique(results_gene2dis_data$gene_symbol[90:180]), mart = ensembl)
chr3_genes <- getBM(attributes=c('hgnc_symbol','chromosome_name','gene_exon_intron','band'), filters =
                       'hgnc_symbol', values =unique(results_gene2dis_data$gene_symbol[180:237]), mart = ensembl)
chr_seq_genes <- rbind(chr1_genes,chr2_genes,chr3_genes)
View(chr_seq_genes)


library(multiMiR)
