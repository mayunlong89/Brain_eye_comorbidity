#########-----------------------------------------
library(ggplot2)
ldsc_original_gwas <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\01-ldsc\\ldsc_original_data.txt",header = TRUE)
ldsc_mtag_gwas <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\01-ldsc\\ldsc_mtag_data.txt",header = TRUE)

##------------Supplementary Figure S12A
ldsc_original_gwas$z_original_GWAS
ldsc_mtag_gwas$z_mtag_GWAS

wilcox.test(ldsc_mtag_gwas$z_mtag_GWAS,ldsc_original_gwas$z_original_GWAS)


ldsc_original_gwas$group <- "original_gwas"
colnames(ldsc_original_gwas) <- c("z_score","group")
ldsc_mtag_gwas$group <- "mtag_gwas"
colnames(ldsc_mtag_gwas) <- c("z_score","group")

data <- rbind(ldsc_original_gwas,ldsc_mtag_gwas)

p<-ggplot(data, aes(x=group, y=z_score,fill=group)) + 
  geom_boxplot()

p+scale_y_continuous(breaks = c(0,10,20,100,300))+
 theme_classic()
