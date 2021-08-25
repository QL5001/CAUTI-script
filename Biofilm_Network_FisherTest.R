#Fisher's exact test of 46 high-biofilm associated genes to prepare data for next-step network analysis using Gehpi software
#Same method applied to the Fisher's exact test of 32 E. coli virulence factors
#Begin#

Bio= read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)  #open the file prepared for Fisher's exact test including the 46 high-biofilm associated genes
BioX= Bio[, 1:2]

BioY= data.frame()
for (i in (1:nrow(BioX)))
{
  BioS= c()
  
  a = c(unlist(BioX[i, ]))
  
  for (j in (1:nrow(BioX)))
  {
    b = c(unlist(BioX[j, ]))
    c = append(a, b)
    d = matrix(c, nrow=2)
    e = fisher.test(d)  
    f = e$p.value
    BioS= c(BioS, f)
    
  }
  
  BioY= rbind(BioY, BioS)
}

write.csv(BioY, "Biofilm_Network_FisherTest.csv")

#Fisher's exact test and pearson correlation coeffcient of 46 high-biofilm associated genes generate data for next-step network analysis using Gephi software
#End#
