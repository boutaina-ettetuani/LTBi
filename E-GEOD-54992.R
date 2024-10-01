
library(affy)
library(oligo)
library(xps)
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(simpleaffy)
library(ArrayExpress)
files_54992 = list.files("E-GEOD-54992", full.names = TRUE)

affy.data_54992 = affy::ReadAffy( filenames = files_54992)
cdfName(affy.data_54992)
####
ph_54992 = affy.data_54992@phenoData
ph_54992
ph_54992$sample
feat_54992 = affy.data_54992@featureData
feat_54992
cdfName(affy.data_54992)
#How to retrieve the IDs of the probe sets that are represented on the arrays
featureNames(affy.data_54992)
#How to retrieve the number of probe sets represented on the arrays 
length(featureNames(affy.data_54992))
#How to retrieve the number of probes represented on the arrays
length(probeNames(affy.data_54992))


#########################################Quality control of microarray data
#Giving the samples informative names
affy.data_54992@phenoData@data
ph_54992[ ,1]

#where we have a data set consisting of two groups of 5/6 replicates, 6 wild type f and 5 M
affy.data_54992@phenoData@data[ ,1] = c("TB1","TB6","TB9","TB5","TB7","TB3","TB13","TB2","TB63","TB8",
                                        "TB93","TB53","TB73","LTBI6","TB4","HD6","HD5","LTBI2","HD1",
                                        "LTBI3","HD3","LTBI4","LTBI1","HD2","HD4","LTBI5","TB16",
                                        "TB96","TB33","TB23","TB66","TB83","TB56","TB43","TB76","TB36","TB26",
                                        "TB86","TB46")


affy.data_54992@phenoData@data

pmexp_54992 = pm(affy.data_54992)

sampleNames_54992 = vector()
logs_54992 = vector()
for (i in 1:11)
{
  sampleNames_54992 = c(sampleNames_54992,rep(affy.data_54992@phenoData@data[i,1],dim(pmexp_54992)[1]))
  logs_54992 = c(logs_54992,log2(pmexp_54992[,i]))
}
#After we have filled the vectors with the data we need, we combine sample names and log intensities into a single data frame:
logData_54992 = data.frame(logInt=logs_54992,sampleName=sampleNames_54992)


###################Calculate quality measures
#########How to calculate the quality measures of each array don't work
data.qc_54992 = qc(affy.data_54992)
avbg(data.qc_54992)
sfs(data.qc_54992)
percent.present(data.qc_54992)
ratios(data.qc)

#####################Normalization of microarray data
#How to normalize the data using RMA 
data.rma_54992 = affy::rma(affy.data_54992)
data.matrix_54992 = exprs(data.rma_54992)

###################################Identification of DE genes
#How to use limma for comparing two groups of samples ?
affy.data_54992@phenoData@data[ ,1] = c("TB","TB","TB","TB","TB","TB","TB","TB","TB","TB",
                                        "TB","TB","TB","LTBI","TB","HD","HD","LTBI","HD",
                                        "LTBI","HD","LTBI","LTBI","HD","HD","LTBI","TB",
                                        "TB","TB","TB","TB","TB","TB","TB","TB","TB","TB",
                                        "TB","TB")




colnames(affy.data_54992@phenoData@data)[1]="source"
affy.data_54992@phenoData@data
groups_54992 = affy.data_54992@phenoData@data$source
f_54992 = factor(groups_54992,levels=c("TB","LTBI","HD"))
design_54992 = model.matrix(~ 0 + f_54992)
colnames(design_54992) = c("TB","LTBI","HD")
data.fit_54992 = lmFit(data.matrix_54992,design_54992)
data.fit_54992$coefficients[1:10,]
contrast.matrix_54992 = makeContrasts(TB-LTBI,TB-HD,LTBI-HD,levels=design_54992)

library(variancePartition)
plotContrasts( contrast.matrix_54992 )

data.fit.con_54992 = contrasts.fit(data.fit_54992,contrast.matrix_54992)
data.fit.eb_54992 = eBayes(data.fit.con_54992)
names(data.fit.eb_54992)
data.fit.eb_54992$coefficients[1:10,]
data.fit.eb_54992$t[1:10,]
data.fit.eb_54992$p.value[1:10,]
topTab_54992<- topTable(data.fit.eb_54992,number = Inf)

topTab_54992

dim(topTab_54992)

#################################
hist(topTab_54992$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Contrasts", xlab = "p-values")


###############################################
#########################################Adjusting for multiple testing and defining DE genes
#How to adjust for multiple testing for a single comparison
options(digits=2)
topTab_54992 = topTable(data.fit.eb,coef=1,number=Inf,adjust.method="BH")



library(hgu133plus2.db)
library(hgu133a.db)
library(hugene10stv1cdf)
library(pd.hugene.2.1.st)



gene_54992_PrbId<- AnnotationDbi::select(hgu133plus2.db,
                                            keys = row.names(gene_54992_0.05),
                                            columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                            keytype = "PROBEID")



length(unique(gene_54992_PrbId$ENTREZID))

anno_subset_54992 <- subset(gene_54992_PrbId, !is.na(SYMBOL))

