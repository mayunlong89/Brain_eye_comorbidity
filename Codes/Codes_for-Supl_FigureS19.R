
###-----------------------------Supplementary Figure S19

mtag_gtex_log_p <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\03-MTAG_LDSC-SEG\\mtag_gtex_logp.txt", header=TRUE)

library(pheatmap)
#COSR_heatmap <- main_results_SCZ$regulon_specificity_s
COSR_heatmap <- mtag_gtex_log_p
COSR_heatmap1 <- COSR_heatmap[,-1]
rownames(COSR_heatmap1) <- COSR_heatmap[,1]

pheatmap::pheatmap(COSR_heatmap1, 
                   cluster_cols = T, 
                   cluster_rows = T,
                   fontsize_row=7.5,
                   scale="column",  # "row", "column" and "none"
                   display_numbers = F) #是否显示每个色块对应的数值(经归一化后的数值)


mat <- cor(COSR_heatmap1)

pheatmap::pheatmap(mat, 
                   cluster_cols = T, 
                   cluster_rows = T,
                   fontsize_row=7.5,
                   #scale="row",  # "row", "column" and "none"
                   display_numbers = T) #是否显示每个色块对应的数值(经归一化后的数值)




mtag_frank_log_p <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\03-MTAG_LDSC-SEG\\mtag_frank_logp.txt", header=TRUE)

library(pheatmap)
#COSR_heatmap <- main_results_SCZ$regulon_specificity_s
mtag_frank_log_p1 <- mtag_frank_log_p
mtag_frank_log_p2 <- mtag_frank_log_p1[,-1]
rownames(mtag_frank_log_p2) <- mtag_frank_log_p1[,1]

pheatmap::pheatmap(mtag_frank_log_p2, 
                   cluster_cols = T, 
                   cluster_rows = T,
                   fontsize_row=3.5,
                   scale="column",  # "row", "column" and "none"
                   display_numbers = F) #是否显示每个色块对应的数值(经归一化后的数值)


pheatmap::pheatmap(mtag_frank_log_p2, 
                   cluster_cols = T, 
                   cluster_rows = T,
                   fontsize_row=4.5,
                   scale="row",  # "row", "column" and "none"
                   display_numbers = F) #是否显示每个色块对应的数值(经归一化后的数值)


mat2 <- cor(mtag_frank_log_p2)

pheatmap::pheatmap(mat2, 
                   cluster_cols = T, 
                   cluster_rows = T,
                   fontsize_row=7.5,
                   #scale="row",  # "row", "column" and "none"
                   display_numbers = T) #是否显示每个色块对应的数值(经归一化后的数值)

