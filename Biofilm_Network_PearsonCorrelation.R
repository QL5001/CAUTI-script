#Pearson correlation coeffcient calculation of 46 high-biofilm associated genes to prepare data for next-step network analysis using Gehpi software
#Same method applied to the Pearson correlation coeffcient calculation of 32 E.coli virulence factors
#Begin#

Bio= read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)  #open the file prepared for Pearson correlation coeffcient calculation including the 46 high-biofilm associated genes
BioX= Bio[, 1:41]

BioY= data.frame()

for (i in (1:nrow(BioX)))
{
  BioS= c()
  
  a = c(unlist(BioX[i, ]))
  
  for (j in (1:nrow(BioX)))
  {
    b = c(unlist(BioX[j, ]))
    c = cor(a, b, method = "pearson")
    BioS= c(BioS, c)
    
  }
  
  BioY= rbind(BioY, BioS)
}

write.csv(BIoY, "Biofilm_Network_-Pearson-result.csv")

#Pearson correlation coeffcient and Fisher's exact test of 46 high-biofilm associated genes generate data for next-step network analysis using Gephi software
#End#
