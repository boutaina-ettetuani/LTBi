
setwd("C:/Users/hp/Desktop/ltbiTranscrp")
library(affy)
library(oligo)
library(xps)
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(simpleaffy)
library(ArrayExpress)
#rawset_41055 = ArrayExpress("E-GEOD-41055")

#rawset_41055 = ArrayExpress("E-GEOD-41055")
#https://wiki.bits.vib.be/index.php/Analyze_your_own_microarray_data_in_R/Bioconductor
files_41055 = list.files("C:/Users/hp/Desktop/ltbiTranscrp/E-GEOD-41055", full.names = TRUE)
################################################How to open CEL files using oligo 
#data_41055 = read.celfiles(files_41055)
#How to open CEL files using affy
affy.data_41055 = read.celfiles( filenames = files_41055)
#bg.corr_41055 <- expresso(affy.data_41055, bg.correct=TRUE, bgcorrect.method="rma",normalize=FALSE, pmcorrect.method="pmonly",summary.method="avgdiff")

#How to retrieve the name of the CDF file associated with the arrays 
cdfName(affy.data_41055)
####
ph_41055 = affy.data_41055@phenoData
ph_41055
ph_41055$sample
feat_41055 = affy.data_41055@featureData
feat_41055
cdfName(affy.data_41055)
#How to retrieve the IDs of the probe sets that are represented on the arrays
featureNames(affy.data_41055)
#How to retrieve the number of probe sets represented on the arrays 
length(featureNames(affy.data_41055))
#How to retrieve the number of probes represented on the arrays
length(probeNames(affy.data_41055))


#########################################Quality control of microarray data
#Giving the samples informative names
affy.data_41055@phenoData@data
ph_41055[ ,1]

#where we have a data set consisting of two groups of 5/6 replicates, 6 wild type f and 5 M
affy.data_41055@phenoData@data[ ,1] = c("healthycontrol1","healthycontrol2","healthycontrol3","healthycontrol4",
                                        "healthycontrol5","healthycontrol6","healthycontrol7","healthycontrol8",
                                        "healthycontrol9",
                                        "latentTBinfection1","latentTBinfection2","latentTBinfection3","latentTBinfection4",
                                        "latentTBinfection5","latentTBinfection6","latentTBinfection7","latentTBinfection8",
                                        "latentTBinfection9",
                                        "activeTBinfection1","activeTBinfection2","activeTBinfection3","activeTBinfection4",
                                        "activeTBinfection5","activeTBinfection6","activeTBinfection7","activeTBinfection8",
                                        "activeTBinfection9")

affy.data_41055@phenoData@data
#How to print the raw intensities of each array
#for (i in 1:11)
#{
#  name = paste("image",i,".jpg",sep="")
#  jpeg(name)
#  image(affy.data_41055[,i],main=affy.data_41055@phenoData@data$sample[i])
#  dev.off()
#}

#op = par(mfrow = c(4,4))
#for (i in 1:11){image(affy.data_41055[,i],main=affy.data_41055@phenoData@data$sample[i])}

#How to create a plot containing the histograms of all samples (each sample in a different color) using ggplot ?
pmexp_41055 = pm(affy.data_41055)

sampleNames_41055 = vector()
logs_41055 = vector()
for (i in 1:11)
{
  sampleNames_41055 = c(sampleNames_41055,rep(affy.data_41055@phenoData@data[i,1],dim(pmexp_41055)[1]))
  logs_41055 = c(logs_41055,log2(pmexp_41055[,i]))
}
#After we have filled the vectors with the data we need, we combine sample names and log intensities into a single data frame:
logData_41055 = data.frame(logInt=logs_41055,sampleName=sampleNames_41055)
#Now we can create the plot:
dataHist2_41055 = ggplot(logData_41055, aes(logInt, colour = sampleName)) 
dataHist2_41055 + geom_density()
##How to create boxplots of microarray data
#How to plot the box plots of the raw data of each array ?
name = "boxplot.jpg"
jpeg(name)
boxplot(affy.data_41055,which='pm',col='red',names=affy.data_41055@phenoData@data$sample) 
dev.off()
#How to create a boxplot of normalized intensities ?
pmexp_41055 = pm(affy.data_41055)
sampleNames_41055 = vector()
logs_41055 = vector()
for (i in 1:11)
{
  sampleNames_41055 = c(sampleNames,rep(affy.data_41055@phenoData@data[i,1],dim(pmexp_41055)[1]))
  logs_41055 = c(logs_41055,log2(pmexp_41055[,i]))
}

logData_41055 = data.frame(logInt=logs_41055,sampleName=sampleNames_41055)
dataBox_41055 = ggplot(logData_41055,aes(sampleName,logInt))
dataBox_41055 + geom_boxplot()
###################Calculate quality measures
#########How to calculate the quality measures of each array don't work
data.qc_41055 = qc(affy.data_41055)
avbg(data.qc_41055)
sfs(data.qc_41055)
percent.present(data.qc_41055)
ratios(data.qc)

#####################Normalization of microarray data
#How to normalize the data using RMA 
data.rma_41055 = rma(affy.data_41055)
data.matrix_41055 = exprs(data.rma_41055)
#####################Checking the effect of the normalization
###How to plot the MA plot of each array
for (i in 1:11)
{
  name = paste("MAplot",i,".jpg",sep="")
  jpeg(name)
  MAplot(affy.data_41055,which=i)
  dev.off()
}
##How to create an MA plot of the normalized intensities
for (i in 1:11)
{
  name = paste("MAplotnorm",i,".jpg",sep="")
  jpeg(name)
  MAplot(data.rma_41055,which=i)
  dev.off()
}

##PCA plot
color=c('green','green','green','green','green','green','blue','blue','blue','blue','blue')
data.PC = prcomp(t(data.matrix_41055),scale.=TRUE)
plot(data.PC$x[1:2],col=color)

###################################Identification of DE genes
#How to use limma for comparing two groups of samples ?
#affy.data_41055@phenoData@data[ ,1] = c("EXLV","EXLV","EXLV","TB","TB","TB","TB","TB",
                                        #"TB","EXLV","EXLV","EXLV","TB","TB","TB","TB","TB",
                                        #"TB","EXLV","EXLV","EXLV","TB","TB","TB","TB","TB","TB")
  
  
affy.data_41055@phenoData@data[ ,1] = c("healthycontrol","healthycontrol","healthycontrol","healthycontrol", "healthycontrol","healthycontrol","healthycontrol","healthycontrol",
                                        "healthycontrol",
                                        "latentTBinfection","latentTBinfection","latentTBinfection","latentTBinfection",
                                       "latentTBinfection","latentTBinfection","latentTBinfection","latentTBinfection",
                                        "latentTBinfection",
                                        "activeTBinfection","activeTBinfection","activeTBinfection","activeTBinfection",
                                        "activeTBinfection","activeTBinfection","activeTBinfection","activeTBinfection",
                                        "activeTBinfection")
colnames(affy.data_41055@phenoData@data)[1]="source"
affy.data_41055@phenoData@data
groups_41055 = affy.data_41055@phenoData@data$source
f_41055 = factor(groups_41055,levels=c("healthycontrol","latentTBinfection","activeTBinfection"))
#f_41055 = factor(groups_41055,levels=c("EXLV","TB"))
design_41055 = model.matrix(~ 0 + f_41055)
colnames(design_41055) = c("healthycontrol","latentTBinfection","activeTBinfection")
data.fit_41055 = lmFit(data.matrix_41055,design_41055)
data.fit_41055$coefficients[1:10,]
contrast.matrix_41055 = makeContrasts(healthycontrol-latentTBinfection,healthycontrol-activeTBinfection,latentTBinfection-activeTBinfection,levels=design_41055)

library(variancePartition)
plotContrasts( contrast.matrix_41055 )

data.fit.con_41055 = contrasts.fit(data.fit_41055,contrast.matrix_41055)
data.fit.eb_41055 = eBayes(data.fit.con_41055)
names(data.fit.eb_41055)
data.fit.eb_41055$coefficients[1:10,]
data.fit.eb_41055$t[1:10,]
data.fit.eb_41055$p.value[1:10,]
topTab_41055<- topTable(data.fit.eb_41055,number = Inf)
#topTab_41055 <-topTable(data.fit.eb_41055,coef=2,number=200,adjust.method="BH")
topTab_41055


#name = "Volcano.jpg"
#jpeg(name)
#volcanoplot(data.fit.eb,coef=1,highlight=10)
#dev.off()
#################################
hist(topTab_41055$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Contrasts", xlab = "p-values")

volcano_names_41055 <- ifelse(abs(data.fit.eb_41055$p.value)<= 0.4, data.fit.eb_41055$coefficients, NA)
#######################################
volcanoplot(data.fit.eb_41055, coef = 1L, style = "p-value", highlight = 100, 
            names = volcano_names_41055, 
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

volcanoplot(data.fit.eb_41055,coef=1,highlight=60,names=data.fit.eb_41055$coefficients,  col="black")
volcanoplot(data.fit.eb_41055,coef=1,highlight=100,  col="black")
###################################################
######################################################################
nrow(subset(topTab_41055, P.Value <= 0.000001))
nrow(subset(topTab_41055, P.Value <= 0.00001))
nrow(subset(topTab_41055, P.Value <= 0.0001))
nrow(subset(topTab_41055, P.Value <= 0.001))
###############################################
gene_41055_0.000001 = subset(topTab_41055, P.Value <=  0.000001)
unique(gene_41055_0.000001)
gene_41055_0.00001 = subset(topTab_41055, P.Value <=  0.00001)
unique(gene_41055_0.00001)
gene_41055_0.0001 = subset(topTab_41055, P.Value <=  0.0001)
unique(gene_41055_0.0001)
gene_41055_0.001 = subset(topTab_41055, P.Value <=  0.001)
unique(gene_41055_0.001)

gene_41055_0.005 = subset(topTab_41055, P.Value <=  0.05)
unique(gene_41055_0.005)
###############################################
#########################################Adjusting for multiple testing and defining DE genes
#How to adjust for multiple testing for a single comparison
options(digits=2)
topTab_41055 = topTable(data.fit.eb,coef=1,number=Inf,adjust.method="BH")

affy.data_41055


library(hgu133plus2.db)
library(hgu133a.db)
library(hugene10stv1cdf)
library(pd.hugene.2.1.st)

library(pd.huex.1.0.st.v2)
library(hugene21sttranscriptcluster.db)




library(biomaRt)
listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl", mart)

gene_41055_pre <- getBM(attributes=c("affy_hg_u133a"),
                    filters = 'affy_huex_1_0_st_v2',
                    values = row.names(gene_41055_0.005),  
                    mart = mart,
                    useCache = FALSE)


gene_41055_PrbId <- AnnotationDbi::select(hgu133plus2.db,
                                     keys = gene_41055_pre$affy_hg_u133a,
                                     columns = c("SYMBOL", "GENENAME","ENSEMBL","ENTREZID"),
                                     keytype = "PROBEID")
gene_41055_PrbId
length(unique(gene_41055_PrbId$ENTREZID))
anno_subset_41055 <- subset(gene_41055_PrbId, !is.na(SYMBOL))





