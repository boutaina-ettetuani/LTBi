
library(affy)
library(oligo)
library(xps)
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(simpleaffy)
library(ArrayExpress)
library(limma)
files_59184 = list.files("E-GEOD-27984", full.names = TRUE)


affy.data_59184 = affy::ReadAffy( filenames = files_59184)

cdfName(affy.data_59184)
####
ph_59184 = affy.data_59184@phenoData
ph_59184
ph_59184$sample
feat_59184 = affy.data_59184@featureData
feat_59184
cdfName(affy.data_59184)
#How to retrieve the IDs of the probe sets that are represented on the arrays
featureNames(affy.data_59184)
#How to retrieve the number of probe sets represented on the arrays 
length(featureNames(affy.data_59184))
#How to retrieve the number of probes represented on the arrays
length(probeNames(affy.data_59184))


#########################################Quality control of microarray data
#Giving the samples informative names
affy.data_59184@phenoData@data
ph_59184[ ,1]

#where we have a data set consisting of two groups of 5/6 replicates, 6 wild type f and 5 M
affy.data_59184@phenoData@data[ ,1] = c("A2","A3","A4","A5","A6","A7","A8","A9",
                                        "B2","B3","B4","B5","B6","B7","B8","B9",
                                        "C2","C3","C4","C5","C6","C7","C8","C9",
                                        "D2","D3","D4","D5","D6","D7","D8","D9")
affy.data_59184@phenoData@data


#How to create a plot containing the histograms of all samples (each sample in a different color) using ggplot ?
pmexp_59184 = pm(affy.data_59184)

###################Calculate quality measures
#########How to calculate the quality measures of each array don't work
data.qc_59184 = qc(affy.data_59184)
avbg(data.qc_59184)
sfs(data.qc_59184)
percent.present(data.qc_59184)
ratios(data.qc)

#####################Normalization of microarray data
#How to normalize the data using RMA 
data.rma_59184 = affy::rma(affy.data_59184)
data.matrix_59184 = exprs(data.rma_59184)


###################################Identification of DE genes
#How to use limma for comparing two groups of samples ?
affy.data_59184@phenoData@data[ ,1] =  c("A","A","A","A","A","A","A","A",
                                         "B","B","B","B","B","B","B","B",
                                         "C","C","C","C","C","C","C","C",
                                         "D","D","D","D","D","D","D","D")
  
colnames(affy.data_59184@phenoData@data)[1]="source"
affy.data_59184@phenoData@data
groups_59184 = affy.data_59184@phenoData@data$source
f_59184 = factor(groups_59184,levels=c("A","B","C","D"))
design_59184 = model.matrix(~ 0 + f_59184)
colnames(design_59184) = c("A","B","C","D")
data.fit_59184 = lmFit(data.matrix_59184,design_59184)
data.fit_59184$coefficients[1:10,]
contrast.matrix_59184 = makeContrasts(A-B,A-C,A-D,B-C,B-D,C-D,levels=design_59184)

library(variancePartition)
plotContrasts( contrast.matrix_59184 )

data.fit.con_59184 = contrasts.fit(data.fit_59184,contrast.matrix_59184)
data.fit.eb_59184 = eBayes(data.fit.con_59184)
names(data.fit.eb_59184)
data.fit.eb_59184$coefficients[1:10,]
data.fit.eb_59184$t[1:10,]
data.fit.eb_59184$p.value[1:10,]
topTab_59184<- topTable(data.fit.eb_59184,number = Inf)

topTab_59184


#################################
hist(topTab_59184$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Contrasts", xlab = "p-values")

volcano_names_59184 <- ifelse(abs(data.fit.eb_59184$p.value)<= 0.4, data.fit.eb_59184$coefficients, NA)
#######################################

###############################################
#########################################Adjusting for multiple testing and defining DE genes
#How to adjust for multiple testing for a single comparison
options(digits=2)
#topTab_59184 = topTable(data.fit.eb,coef=1,number=Inf,adjust.method="BH")


library(hgu133plus2.db)
library(hgu133a.db)
library(hugene10stv1cdf)
library(pd.hugene.2.1.st)



gene_59184_PrbId<- AnnotationDbi::select(hgu133plus2.db,
                                         keys = row.names(gene_59184_0.05),
                                         columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                         keytype = "PROBEID")



length(unique(gene_59184_PrbId$ENTREZID))

anno_subset_59184 <- subset(gene_59184_PrbId, !is.na(SYMBOL))
