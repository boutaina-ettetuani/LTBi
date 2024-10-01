
union(setA,setB)

lit = c("FCGR1A","HK3","RAB13","RBBP8","IFI44L","TIMM10","BCL6",
        "SMARCD3","CYP4F3","SLPI","PTPRC","ASUN","DHX29","NEFM",
        "IGF2","IGFALS","IGFBP3","DEFA1","GBP5","GBP6","CHRM2",
        "SNX17","AMPH","PIGC","TAS2R46","HBD","GLDC","ACOT7",
        "S100P","STYXL1","FRP1","DEFA1","DEFA3","DEFA4","BPI",
        "CD64","LTF","Rab33A","STAT1","IFIT2","GBP5","STAT2",
        "GBP3","GBP2","NMI","IFIT2","GSDMD","ORL1","TCN1","TCN2")

lit_entrezid <-  AnnotationDbi::select(hgu133plus2.db,
                                       keys = lit,
                                       columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                        keytype = "SYMBOL")



length(unique(gene_59184_PrbId$ENTREZID))
EGEOD54992= unique(anno_subset_54992$SYMBOL)
EGEOD41055 = unique(anno_subset_41055$SYMBOL)
EGEOD27984 =unique(anno_subset_59184$SYMBOL)
intersect(lit,EGEOD54992)


tst <- c(unique(lit),unique(EGEOD54992),unique(EGEOD41055),unique(EGEOD27984))
tst <- tst[duplicated(tst)]
tst[duplicated(tst)]

library("dplyr")  
datafull_join <- rbind(anno_subset_54992[,c("SYMBOL","ENTREZID")], anno_subset_41055[,c("SYMBOL","ENTREZID")],anno_subset_59184[,c("SYMBOL","ENTREZID")],lit_entrezid[,c("SYMBOL","ENTREZID")])


library(systemPipeR)



setlist_symbol <- list(lit= lit, EGEOD54992= unique(anno_subset_54992$SYMBOL),EGEOD41055 = unique(anno_subset_41055$SYMBOL), EGEOD27984 =unique(anno_subset_59184$SYMBOL))
## 2-way Venn diagram
vennset2_symbol <- overLapper(setlist_symbol[1:4], type="vennsets")
#vennset2_symbol <- overLapper(setlist_symbol[c(1,3)], type="vennsets")

vennPlot(vennset2_symbol)


library(VennDiagram)
library(purrr)
library(RVenn)
library(ggplot2)
toy <- list(lit= lit,EGEOD54992= unique(anno_subset_54992$SYMBOL),EGEOD41055 = unique(anno_subset_41055$SYMBOL), EGEOD27984 =unique(anno_subset_59184$SYMBOL))
toy = Venn(toy)
overlap(toy)
unite(toy)
overlap(toy, c("lit", "EGEOD54992", "EGEOD41055", "EGEOD27984"))
overlap_pairs(toy, slice = 1:4)
#####Draw the Venn Diagram
ggvenn(toy, slice = c(1, 2,4))

ggvenn(toy, slice = c(1, 2, 3,4),fill = c("deeppink", "dodgerblue3", "pink"))
#fill = c("gold", "dodgerblue3", "deeppink")
ggvenn(toy, slice = c(1, 2, 3,4),fill = c("deeppink", "dodgerblue3", "pink"))
