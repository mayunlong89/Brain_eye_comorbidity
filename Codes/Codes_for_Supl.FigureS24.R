#############--------------------------------Heatmap for the fourth pleitropic loci----cis-eQTL
data_fourth_loci <- read.csv("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\04-metasoft\\04-森林图---PM-plot\\Top-4\\PACSIN2_cisQTL_data.csv")
library(pheatmap)
data_fourth_loci1 <- data_fourth_loci
data_fourth_loci2 <- data_fourth_loci1[,-1]
rownames(data_fourth_loci2) <- data_fourth_loci1[,1]
pheatmap::pheatmap(data_fourth_loci2, 
                   cluster_cols = T, 
                   cluster_rows = T,
                   fontsize_row=7.5,
                   #scale="row",  # "row", "column" and "none"
                   display_numbers = F) #是否显示每个色块对应的数值(经归一化后的数值)
