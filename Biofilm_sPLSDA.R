#Sparse partial least squares discriminant analysis (sPLSDA) of the whole genome sequences of 29 E. coli isolates including 13 high and 16 low biofilm formers
#Begin#

library(knitr)
library(mixOmics)
library(snow)

WGS = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE) #open the whole genome seqeunce (WGS) of 29 E. coli isolates including 13 high and 16 low biofilm formers

X = WGS[, 1:15993]
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y = WGS$Biofilm

splsda.WGS = splsda(X, Y, ncomp = 10)
barplot(splsda.WGS$explained_variance$X, ylim = c(0, 0.15), las=2, cex.names=1,  main = "Biofilm_sPLSDA_ExplainedVarianceHistogram")

ev = splsda.WGS$explained_variance   #explained_variance
write.csv(ev, "Biofilm_sPLSDA_ExplainedVariance.csv")

sco = splsda.WGS$variates   #prinicipal components
write.csv(sco, "Biofilm_sPLSDA_PrincipleComponent.csv")

loa = data.frame(splsda.WGS$loadings$X)   #loading values
gene_name = colnames(X)   #gene names

WGS_1 = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)   #open the whole genome seqeunce (WGS) of 29 low and high biofilm E. coli isolates to extract annotated gene functions
gene_anno = WGS_1$Annotation   #gene functions

loa_PC1 = data.frame(gene = gene_name, annotation = gene_anno, PC1loa = loa[, 1])
write.csv(loa_PC1, "Biofilm_sPLSDA_loading_PC1.csv")

loa_PC2 = data.frame(gene = gene_name, annotatoin = gene_anno, PC2loa = loa[, 2])
write.csv(loa_PC2, "Biofilm_sPLSDA_loading_PC2.csv")

par(mar=c(5,4,4,2)+0.1)
png(filename = "Biofilm_sPLSDA_ScorePlot.png", width = 600, height = 400)  #sPLSDA score plot
plotIndiv(splsda.WGS, comp = c(1,2), group = WGS$Biofilm, pch = c(15, 16), cex = 3, col = rainbow(2), ellipse = FALSE, ind.names = FALSE, legend = TRUE, size.legend = 15, size.xlabel = 18, size.ylabel = 18, size.axis = 15, X.label = "Component 1(9%)", Y.label = "Component 2(10%)", title = "Biofilm_sPLSDA_ScorePlot")
graphics.off()

#Exported csv files are used for further data analyses and drawing figures in Excel or Prism
#End#
