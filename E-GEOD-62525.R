library(affy)
library(oligo)
library(xps)
library(Matrix)
library(lattice)
library(fdrtool)
library(rpart)
library(simpleaffy)
library(ArrayExpress)
rawset_62525 = ArrayExpress("E-GEOD-62525")


files_62525 = list.files("E-GEOD-62525", full.names = TRUE)
################################################How to open CEL files using oligo 

affy.data_62525 = affy::ReadAffy( filenames = files_62525)

#How to retrieve the name of the CDF file associated with the arrays 
cdfName(affy.data_62525)
####
ph_62525 = affy.data_62525@phenoData
ph_62525
ph_62525$sample
feat_62525 = affy.data_62525@featureData
feat_62525
cdfName(affy.data_62525)
#How to retrieve the IDs of the probe sets that are represented on the arrays
featureNames(affy.data_62525)
#How to retrieve the number of probe sets represented on the arrays 
length(featureNames(affy.data_62525))
#How to retrieve the number of probes represented on the arrays
length(probeNames(affy.data_62525))


#########################################Quality control of microarray data
#Giving the samples informative names
affy.data_62525@phenoData@data
ph_62525[ ,1]

#where we have a data set consisting of two groups of 5/6 replicates, 6 wild type f and 5 M
affy.data_62525@phenoData@data[ ,1] = c("glomeruli_FSGS_repF1","glomeruli_FSGS_repF2","glomeruli_FSGS_repF3","glomeruli_FSGS_repF4",
                                        "glomeruli_FSGS_repF5","glomeruli_FSGS_repF6","glomeruli_FSGS_repM1","glomeruli_FSGS_repM2",
                                        "glomeruli_FSGS_repM3","glomeruli_FSGS_repM4","glomeruli_FSGS_repM5")
affy.data_62525@phenoData@data


pmexp_62525 = pm(affy.data_62525)

###################Calculate quality measures
#########How to calculate the quality measures of each array don't work
data.qc_62525 = qc(affy.data_62525)
avbg(data.qc_62525)
sfs(data.qc_62525)
percent.present(data.qc_62525)
ratios(data.qc)

#####################Normalization of microarray data
#How to normalize the data using RMA 
data.rma_62525 = rma(affy.data_62525)
data.matrix_62525 = exprs(data.rma_62525)

###################################Identification of DE genes
#How to use limma for comparing two groups of samples ?
affy.data_62525@phenoData@data[ ,1] = c("repF","repF","repF","repF","repF","repF","repM","repM","repM","repM","repM")
colnames(affy.data_62525@phenoData@data)[1]="source"
affy.data_62525@phenoData@data
groups_62525 = affy.data_62525@phenoData@data$source
f_62525 = factor(groups_62525,levels=c("repF","repM"))
design_62525 = model.matrix(~ 0 + f_62525)
colnames(design_62525) = c("repF","repM")
data.fit_62525 = lmFit(data.matrix_62525,design_62525)
data.fit_62525$coefficients[1:10,]
contrast.matrix_62525 = makeContrasts(repF-repM,levels=design_62525)

library(variancePartition)
plotContrasts( contrast.matrix_62525 )

data.fit.con_62525 = contrasts.fit(data.fit_62525,contrast.matrix_62525)
data.fit.eb_62525 = eBayes(data.fit.con_62525)
names(data.fit.eb_62525)
data.fit.eb_62525$coefficients[1:10,]
data.fit.eb_62525$t[1:10,]
data.fit.eb_62525$p.value[1:10,]
topTab_62525<- topTable(data.fit.eb_62525,number = Inf)
#topTab_62525 <-topTable(data.fit.eb_62525,coef=2,number=200,adjust.method="BH")
topTab_62525

group2 <- group
levels(group2) <- c("basal.lactate","basal.preg","basal.virgin","lum.lactate", "lum.preg", "lum.virgin")

#name = "Volcano.jpg"
#jpeg(name)
#volcanoplot(data.fit.eb,coef=1,highlight=10)
#dev.off()
#################################
hist(topTab_62525$P.Value, col = brewer.pal(3, name = "Set2")[1],
     main = "Contrasts", xlab = "p-values")

volcano_names_62525 <- ifelse(abs(data.fit.eb_62525$p.value)<= 0.4, data.fit.eb_62525$coefficients, NA)
#######################################
volcanoplot(data.fit.eb_62525, coef = 1L, style = "p-value", highlight = 100, 
            names = volcano_names_62525, 
            xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35)

volcanoplot(data.fit.eb_62525,coef=1,highlight=60,names=data.fit.eb_62525$coefficients,  col="black")

###############################################
#########################################Adjusting for multiple testing and defining DE genes
#How to adjust for multiple testing for a single comparison
options(digits=2)
topTab_62525 = topTable(data.fit.eb,coef=1,number=Inf,adjust.method="BH")

