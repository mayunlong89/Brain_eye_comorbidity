#######----------the Correlations across Five factors estimated by Genomic SEM 
data_factor_cor <- read.table("./02-GenomicSEM/factor_cor.txt",header = TRUE)

library(pheatmap)
#COSR_heatmap <- main_results_SCZ$regulon_specificity_s
COSR_heatmap <- data_factor_cor
COSR_heatmap1 <- COSR_heatmap[,-1]
rownames(COSR_heatmap1) <- COSR_heatmap[,1]
pheatmap::pheatmap(COSR_heatmap1, 
                   cluster_cols = T, 
                   cluster_rows = T,
                   fontsize_row=9.5,
                   display_numbers = TRUE) #是否显示每个色块对应的数值(经归一化后的数值)
