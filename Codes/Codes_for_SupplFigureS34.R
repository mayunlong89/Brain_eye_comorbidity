

#####------------------Heatmap for the addmodulescore of single-cell data
## eye-brain pleiotropic genes vs. non-pleiotropic genes


headmap_addm <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\05-singlecell\\heatmap_addmodulescore.txt",header=TRUE)

library(pheatmap)

headmap_addm1 <- headmap_addm[,-1]
rownames(headmap_addm1)<-headmap_addm[,1]
xx <- t(headmap_addm1)
pheatmap::pheatmap(xx, 
                   cluster_cols = F, 
                   cluster_rows = F,
                   fontsize_row=7.5,
                   scale="none",  # "row", "column" and "none"
                   display_numbers = F) #是否显示每个色块对应的数值(经归一化后的数值)


