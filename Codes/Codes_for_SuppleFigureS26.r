
###---------------Brain-based single-cell data
# Load Packages
library(Seurat)

brain_data <- readRDS("F:\\BaiduNetdiskDownload\\眼脑单细胞\\resolution0.2_final4.rds")



###---------------SMultiXcan-MAGMA correlation

data_cor_gene <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\04-metasoft\\09-H-MAGMA-S-MultiXcan\\magma_smultixcan_gene_count_cor.txt", header=T)


cor.test(data_cor_gene$MAGMA,data_cor_gene$SMultiXcan)

plot(data_cor_gene$MAGMA~data_cor_gene$SMultiXcan)
abline(data_cor_gene$MAGMA~data_cor_gene$SMultiXcan)

library(ggplot2) ##s-------------------------------------------Scatter plot-smooth line
#library(ggpmisc) #install.packages("ggpmisc")
ggplot(data=data_cor_gene,mapping = aes(x=MAGMA,y=SMultiXcan))+
  geom_point(pch=19,col= "#69B9DA",alpha=1,cex=3.0)+  #alpha 设置透明度，越小越透明
  geom_smooth(method = "lm",formula = y~x, col="red",alpha=0.2)+
  theme_classic()  ###这个放到前面不要放到最后，否则转角angle不起作用



###---------------SMultiXcan-MAGMA correlation

data_cor_disease <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\04-metasoft\\09-H-MAGMA-S-MultiXcan\\magma_smultixcan_disease_count_cor.txt", header=T)


cor.test(data_cor_disease$Total_smultixcan,data_cor_disease$Total_magma)
cor.test(data_cor_disease$Porportion_smultixcan,data_cor_disease$Porportion_magma)

library(ggplot2) ##s-------------------------------------------Scatter plot-smooth line
#library(ggpmisc) #install.packages("ggpmisc")
ggplot(data=data_cor_disease,mapping = aes(x=Porportion_smultixcan,y=Porportion_magma))+
  geom_point(pch=19,col= "#69B9DA",alpha=1,cex=3.0)+  #alpha 设置透明度，越小越透明
  geom_smooth(method = "lm",formula = y~x, col="red",alpha=0.2)+
  theme_classic()  ###这个放到前面不要放到最后，否则转角angle不起作用

