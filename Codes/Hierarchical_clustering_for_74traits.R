#Hierarchical clustering 
library(pvclust)  #install.packages("pvclust")
traits_74_cor <- read.csv("./01-LDSC_other74_traits/rg.csv")

traits_74_cor1 <- traits_74_cor[,-1]
rownames(traits_74_cor1) <- traits_74_cor[,1]

result <- pvclust(traits_74_cor1, method.dist = "cor", method.hclust = "ward.D",
                  nboot = 1000, parallel = T)
plot(result)
pvrect(result, alpha = 0.95)
