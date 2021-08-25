#Sparse principal component analysis (sPCA) of the whole genome sequences of phylotype B2 E. coli isolates
#Begin#

library(knitr)
library(mixOmics)
library(dplyr)
library(stats)

WGS = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)    #open the whole genome seqeunce (WGS) of phylotype B2 isolates

X = WGS[, 1:15993]  
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y = WGS$B2Branch

spca.WGS = spca(X_center, ncomp = 10, center = FALSE, scale = FALSE)
plot(pca.WGS, ylim = c(0, 0.5), main = "B2Phylotype_sPCA_ExplainedVarianceHistogram")

ev = spca.WGS$explained_variance   #explained variance
write.csv(ev, "B2Phylotype_sPCA_ExplainedVariance.csv")

sco = spca.WGS$variates   #prinicipal components
write.csv(sco, "B2Phylotype_sPCA__PrincipleComponent.csv")

loa = data.frame(spca.WGS$loadings$X)   #loading values
gene_name = colnames(X)   #gene names
WGS_1 = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)   #open the whole genome seqeunce (WGS) of phylotype B2 isolates to extract annotated gene functions
gene_anno = WGS_1$Annotation   #gene functions

loa_PC1 = data.frame(gene = gene_name, annotation = gene_anno, PC1loa = loa[, 1]) #PC1 loading values
write.csv(loa_PC1, "B2Phylotype_sPCA_loading_PC1.csv")

loa_PC2 = data.frame(gene = gene_name, annotatoin = gene_anno, PC2loa = loa[, 2]) #PC2 loading values
write.csv(loa_PC2, "B2Phylotype_sPCA_loading_PC2.csv")

par(mar=c(5,4,4,2)+0.1)
png(filename = "B2Phylotype_sPCA_ScorePlot.png", width = 600, height = 400)  #sPCA score plot
plotIndiv(spca.WGS, comp = c(1,2), group = WGS$B2Branch, pch = c(15, 16), cex = 4, col = rainbow(2), ellipse = FALSE, ind.names = FALSE, legend = TRUE, size.legend = 15, size.xlabel = 18, size.ylabel = 18, X.label = "PC1(31%)", Y.label = "PC2(10%)", size.axis = 15, title = "B2Phylotype_sPCA_ScorePlot")
graphics.off()

#Exported csv files are used for further data analyses and drawing figures in Excel or Prism
#End#
