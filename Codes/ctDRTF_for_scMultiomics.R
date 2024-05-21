
#@ 2024-05-14

#@ ctDRTF for human brain single-cell multiomic data for eye-brain connections
suppressPackageStartupMessages({
  library(Pando)
  library(Seurat)
  library(Signac)
  library(ArchR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(dplyr)
  library(dotgen)
  library(ctDRTF)
})




set.seed(1234)

pbmc_10x_brain <- readRDS("/Users/mayunlong/Desktop/WMU2024/ctDRTF-codes and data/scHBO/10X_human_brain/10X_Human_Brain.rds")

pbmc_10x_brain$cell_type

table(pbmc_10x_brain$cell_type)

Idents(pbmc_10x_brain) <- pbmc_10x_brain$cell_type

DimPlot(pbmc_10x_brain)



##running for human brain single-cell multiomic data
grn_outputs <- GRN_func(single_cell = pbmc_10x_brain,n_genes=10)


saveRDS(grn_outputs,file = "~/Desktop/human_brain_2464cells_10x_scMultioimics_grn.rds") 


# calculate specificity scores
Idents(pbmc_10x_brain) <- pbmc_10x_brain$cell_type
data_s1 <- COSR_pre_func(single_cell = pbmc_10x_brain)




#Eight psychiatric disorders-----MAGMA results

setwd("~/Desktop/WMU2024/ctDRTF-codes and data/scHBO/GWAS-data/MAGMA_result_fuma/eight_psychiatric_disease_fuma/")

#ADHD
data_ADHD <- read.table("magma.genes_ADHD.out",header = TRUE)
data_ADHD <- data_ADHD %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_ADHD <- data_ADHD[,c(10,11,8)] 
head(magma_GWAS_ADHD)

output_final1 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_ADHD,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final1,file="ctDRTF_v1_ADHD_1000times_scaled_adjusted_mi=1_top500genes.csv")




#AN
data_AN <- read.table("magma.genes_AN.out",header = TRUE)
data_AN <- data_AN %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_AN <- data_AN[,c(10,11,8)] 
head(magma_GWAS_AN)

output_final2 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_AN,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final2,file="ctDRTF_v1_AN_1000times_scaled_adjusted_mi=1_top500genes.csv")





#ASD
data_ASD <- read.table("magma.genes_ASD.out",header = TRUE)
data_ASD <- data_ASD %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_ASD <- data_ASD[,c(10,11,8)] 
head(magma_GWAS_ASD)


output_final3 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_ASD,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final3,file="ctDRTF_v1_ASD_1000times_scaled_adjusted_mi=1_top500genes.csv")




#BIP
data_BIP <- read.table("magma.genes_BIP.out",header = TRUE)
data_BIP <- data_BIP %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_BIP <- data_BIP[,c(10,11,8)] 
head(magma_GWAS_BIP)

output_final4 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_BIP,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final4,file="ctDRTF_v1_BIP_1000times_scaled_adjusted_mi=1_top500genes.csv")


#MDD
data_MDD <- read.table("magma.genes_MDD.out",header = TRUE)
data_MDD <- data_MDD %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_MDD <- data_MDD[,c(10,11,8)] 
head(magma_GWAS_MDD)


output_final5 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_MDD,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final5,file="ctDRTF_v1_MDD_1000times_scaled_adjusted_mi=1_top500genes.csv")


#OCD
data_OCD <- read.table("magma.genes_OCD.out",header = TRUE)
data_OCD <- data_OCD %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_OCD <- data_OCD[,c(10,11,8)] 
head(magma_GWAS_OCD)

output_final6 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_OCD,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final6,file="ctDRTF_v1_OCD_1000times_scaled_adjusted_mi=1_top500genes.csv")

#SCZ
data_SCZ <- read.table("magma.genes_SCZ.out",header = TRUE)
data_SCZ <- data_SCZ %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_SCZ <- data_SCZ[,c(10,11,8)] 
head(magma_GWAS_SCZ)

output_final7 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_SCZ,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final7,file="ctDRTF_v1_SCZ_1000times_scaled_adjusted_mi=1_top500genes.csv")

#TS
data_TS <- read.table("magma.genes_TS.out",header = TRUE)
data_TS <- data_TS %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_TS <- data_TS[,c(10,11,8)] 
head(magma_GWAS_TS)


output_final8 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_TS,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final8,file="ctDRTF_v1_TS_1000times_scaled_adjusted_mi=1_top500genes.csv")

#DR
data_DR <- read.table("magma.genes_DR.out",header = TRUE)
data_DR <- data_DR %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_DR <- data_DR[,c(10,11,8)] 
head(magma_GWAS_DR)


output_final9 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_DR,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final9,file="ctDRTF_v1_DR_1000times_scaled_adjusted_mi=1_top500genes.csv")

#myopia
data_myopia <- read.table("magma.genes_myopia.out",header = TRUE)
data_myopia <- data_myopia %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_myopia <- data_myopia[,c(10,11,8)] 
head(magma_GWAS_myopia)


output_final10 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                  data_regulons1 = grn_outputs$grn,
                                  data_s1 = data_s1,
                                  MAGMA_GWAS_data = magma_GWAS_myopia,
                                  MC_num = 1000,
                                  Gene_num = 500,
                                  theta = 0.5,
                                  mi=1,
                                  mo=1)

write.csv(output_final10,file="ctDRTF_v1_myopia_1000times_scaled_adjusted_mi=1_top500genes.csv")


#AMD
data_AMD <- read.table("magma.genes_AMD.out",header = TRUE)
data_AMD <- data_AMD %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_AMD <- data_AMD[,c(10,11,8)] 
head(magma_GWAS_AMD)

output_final11 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                   data_regulons1 = grn_outputs$grn,
                                   data_s1 = data_s1,
                                   MAGMA_GWAS_data = magma_GWAS_AMD,
                                   MC_num = 1000,
                                   Gene_num = 500,
                                   theta = 0.5,
                                   mi=1,
                                   mo=1)

write.csv(output_final11,file="ctDRTF_v1_AMD_1000times_scaled_adjusted_mi=1_top500genes.csv")


#RD
data_RD <- read.table("magma.genes_RD.out",header = TRUE)
data_RD <- data_RD %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_RD <- data_RD[,c(10,11,8)] 
head(magma_GWAS_RD)

output_final12 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                   data_regulons1 = grn_outputs$grn,
                                   data_s1 = data_s1,
                                   MAGMA_GWAS_data = magma_GWAS_RD,
                                   MC_num = 1000,
                                   Gene_num = 500,
                                   theta = 0.5,
                                   mi=1,
                                   mo=1)

write.csv(output_final12,file="ctDRTF_v1_RD_1000times_scaled_adjusted_mi=1_top500genes.csv")


#GLC
data_GLC <- read.table("magma.genes_GLC.out",header = TRUE)
data_GLC <- data_GLC %>% mutate(logP = -log10(P)) %>% arrange(desc(logP))
magma_GWAS_GLC <- data_GLC[,c(10,11,8)] 
head(magma_GWAS_GLC)


output_final13 <- COSR_func_weight(tf_left = grn_outputs$tf_names,
                                   data_regulons1 = grn_outputs$grn,
                                   data_s1 = data_s1,
                                   MAGMA_GWAS_data = magma_GWAS_GLC,
                                   MC_num = 1000,
                                   Gene_num = 500,
                                   theta = 0.5,
                                   mi=1,
                                   mo=1)

write.csv(output_final13,file="ctDRTF_v1_GLC_1000times_scaled_adjusted_mi=1_top500genes.csv")




##----------Heatmap plots
#@ Heatmap

library(pheatmap)


microglia <- read.table("heatmap_microglia.txt",header = TRUE,sep = "\t")


data_h1 <- microglia[,-1]
rownames(data_h1) <- microglia[,1]


data_h1_score <- data_h1[,c(1,3,5,7,9,11,13,15,17,19,21,23,25)]

data_h1_p <- data_h1[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]

data_h1_p2 <- filter(data_h1_p,ADHD_p < 0.05 | AN_p <0.05 | AMD_p < 0.05 | ASD_p <0.05| DR_p < 0.05 | RD_p < 0.05| BIP_p < 0.05 | GLC_p < 0.05 | MDD_p < 0.05| SCZ_p < 0.05 | TS_p < 0.05 | OCD_p < 0.05| myopia_p < 0.05)
data_h1_score2 <- data_h1_score[which(rownames(data_h1_score) %in% rownames(data_h1_p2)),]

p_vals<-t(data_h1_p2)
heatmap_data <- t(data_h1_score2)

p_vals<-data_h1_p2
heatmap_data <- data_h1_score2

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =4, cellheight =4,fontsize=5,border_color="black",scale="none",
                          fontsize_row=5, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699',"blue")))(102),
                          cluster_rows = T,
                          cluster_cols = T,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data >1.5, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))

##version 2
data_h1_mg2 <- filter(data_h1,(ADHD_p < 0.05 & ADHD_score > 1.5)| (AN_p < 0.05 & AN_score > 1.5) |(AMD_p < 0.05 & AMD_score > 1.5)| (ASD_p < 0.05 & ASD_score > 1.5) | (DR_p < 0.05 & DR_score > 1.5) | (RD_p < 0.05 & RD_score > 1.5) | (BIP_p < 0.05 & BIP_score > 1.5) | (GLC_p < 0.05 & GLC_score  > 1.5) | (MDD_p < 0.05 & MDD_score > 1.5)| (SCZ_p < 0.05 & SCZ_score > 1.5) | (TS_p < 0.05 & TS_score > 1.5) | (OCD_p < 0.05 & OCD_score > 1.5)| (myopia_p < 0.05 & myopia_score > 1.5))

 
data_h1_score2_mg <- data_h1_mg2[,c(1,3,7,9,15,19,23,25,5,11,13,17,21)]

data_h1_p2_mg <- data_h1_mg2[,c(2,4,8,10,16,20,24,26,6,12,14,18,22)]


#p_vals<-t(data_h1_p2_mg)
#heatmap_data <- t(data_h1_score2_mg)

p_vals<-data_h1_p2_mg
heatmap_data <- data_h1_score2_mg

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight =8,fontsize=8,border_color="black",scale="none",
                          fontsize_row=8, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84","white","lightblue",'#336699',"blue")))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data >1.5, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))


 


##exhibitory neurons
ex_neuron <- read.table("heatmap_ex_neurons.txt",header = TRUE,sep = "\t")


data_h1_ex <- ex_neuron[,-1]
rownames(data_h1_ex) <- ex_neuron[,1]


data_h1_ex2 <- filter(data_h1_ex,(ADHD_p < 0.05 & ADHD_score > 1.5)| (AN_p < 0.05 & AN_score > 1.5) |(AMD_p < 0.05 & AMD_score > 1.5)| (ASD_p < 0.05 & ASD_score > 1.5) | (DR_p < 0.05 & DR_score > 1.5) | (RD_p < 0.05 & RD_score > 1.5) | (BIP_p < 0.05 & BIP_score > 1.5) | (GLC_p < 0.05 & GLC_score  > 1.5) | (MDD_p < 0.05 & MDD_score > 1.5)| (SCZ_p < 0.05 & SCZ_score > 1.5) | (TS_p < 0.05 & TS_score > 1.5) | (OCD_p < 0.05 & OCD_score > 1.5)| (myopia_p < 0.05 & myopia_score > 1.5))
 
length(data_h1_ex2$ADHD_score)

#data_h1_score2_ex <- data_h1_ex2[,c(1,3,5,7,9,11,13,15,17,19,21,23,25)]
#data_h1_p2_ex <- data_h1_ex2[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
data_h1_score2_ex <- data_h1_ex2[,c(1,3,7,9,15,19,23,25,5,11,13,17,21)]
data_h1_p2_ex <- data_h1_ex2[,c(2,4,8,10,16,20,24,26,6,12,14,18,22)]


#p_vals<-t(data_h1_p2_ex)
#heatmap_data <- t(data_h1_score2_ex)

p_vals<-data_h1_p2_ex
heatmap_data <- data_h1_score2_ex

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight =8,fontsize=8,border_color="black",scale="none",
                          fontsize_row=8, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77","#E77A77", "#ECBA84","white","lightblue",'#336699',"blue")))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data >1.5, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))





##inhibitory neurons
in_neuron <- read.table("heatmap_in_neurons.txt",header = TRUE,sep = "\t")


data_h1_ex <- in_neuron[,-1]
rownames(data_h1_ex) <- in_neuron[,1]


data_h1_ex2 <- filter(data_h1_ex,(ADHD_p < 0.05 & ADHD_score > 1.5)| (AN_p < 0.05 & AN_score > 1.5) |(AMD_p < 0.05 & AMD_score > 1.5)| (ASD_p < 0.05 & ASD_score > 1.5) | (DR_p < 0.05 & DR_score > 1.5) | (RD_p < 0.05 & RD_score > 1.5) | (BIP_p < 0.05 & BIP_score > 1.5) | (GLC_p < 0.05 & GLC_score  > 1.5) | (MDD_p < 0.05 & MDD_score > 1.5)| (SCZ_p < 0.05 & SCZ_score > 1.5) | (TS_p < 0.05 & TS_score > 1.5) | (OCD_p < 0.05 & OCD_score > 1.5)| (myopia_p < 0.05 & myopia_score > 1.5))

length(data_h1_ex2$ADHD_score)

#data_h1_score2_ex <- data_h1_ex2[,c(1,3,5,7,9,11,13,15,17,19,21,23,25)]
#data_h1_p2_ex <- data_h1_ex2[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
data_h1_score2_ex <- data_h1_ex2[,c(1,3,7,9,15,19,23,25,5,11,13,17,21)]
data_h1_p2_ex <- data_h1_ex2[,c(2,4,8,10,16,20,24,26,6,12,14,18,22)]


#p_vals<-t(data_h1_p2_ex)
#heatmap_data <- t(data_h1_score2_ex)

p_vals<-data_h1_p2_ex
heatmap_data <- data_h1_score2_ex

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight =8,fontsize=8,border_color="black",scale="none",
                          fontsize_row=8, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77","#E77A77", "#ECBA84","white","lightblue",'#336699',"blue")))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data >1.5, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))



##OPC
opc <- read.table("heatmap_opc.txt",header = TRUE,sep = "\t")


data_h1_ex <- opc[,-1]
rownames(data_h1_ex) <- opc[,1]


data_h1_ex2 <- filter(data_h1_ex,(ADHD_p < 0.05 & ADHD_score > 1.5)| (AN_p < 0.05 & AN_score > 1.5) |(AMD_p < 0.05 & AMD_score > 1.5)| (ASD_p < 0.05 & ASD_score > 1.5) | (DR_p < 0.05 & DR_score > 1.5) | (RD_p < 0.05 & RD_score > 1.5) | (BIP_p < 0.05 & BIP_score > 1.5) | (GLC_p < 0.05 & GLC_score  > 1.5) | (MDD_p < 0.05 & MDD_score > 1.5)| (SCZ_p < 0.05 & SCZ_score > 1.5) | (TS_p < 0.05 & TS_score > 1.5) | (OCD_p < 0.05 & OCD_score > 1.5)| (myopia_p < 0.05 & myopia_score > 1.5))

length(data_h1_ex2$ADHD_score)

#data_h1_score2_ex <- data_h1_ex2[,c(1,3,5,7,9,11,13,15,17,19,21,23,25)]
#data_h1_p2_ex <- data_h1_ex2[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
data_h1_score2_ex <- data_h1_ex2[,c(1,3,7,9,15,19,23,25,5,11,13,17,21)]
data_h1_p2_ex <- data_h1_ex2[,c(2,4,8,10,16,20,24,26,6,12,14,18,22)]


#p_vals<-t(data_h1_p2_ex)
#heatmap_data <- t(data_h1_score2_ex)

p_vals<-data_h1_p2_ex
heatmap_data <- data_h1_score2_ex

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight =8,fontsize=8,border_color="black",scale="none",
                          fontsize_row=8, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77","#E77A77", "#ECBA84","#ECBA84","white","lightblue",'#336699',"blue")))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data >1.5, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))



##Astrocyte
astro <- read.table("heatmap_astrocytes.txt",header = TRUE,sep = "\t")


data_h1_ex <- astro[,-1]
rownames(data_h1_ex) <- astro[,1]


data_h1_ex2 <- filter(data_h1_ex,(ADHD_p < 0.05 & ADHD_score > 1.5)| (AN_p < 0.05 & AN_score > 1.5) |(AMD_p < 0.05 & AMD_score > 1.5)| (ASD_p < 0.05 & ASD_score > 1.5) | (DR_p < 0.05 & DR_score > 1.5) | (RD_p < 0.05 & RD_score > 1.5) | (BIP_p < 0.05 & BIP_score > 1.5) | (GLC_p < 0.05 & GLC_score  > 1.5) | (MDD_p < 0.05 & MDD_score > 1.5)| (SCZ_p < 0.05 & SCZ_score > 1.5) | (TS_p < 0.05 & TS_score > 1.5) | (OCD_p < 0.05 & OCD_score > 1.5)| (myopia_p < 0.05 & myopia_score > 1.5))

length(data_h1_ex2$ADHD_score)

#data_h1_score2_ex <- data_h1_ex2[,c(1,3,5,7,9,11,13,15,17,19,21,23,25)]
#data_h1_p2_ex <- data_h1_ex2[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
data_h1_score2_ex <- data_h1_ex2[,c(1,3,7,9,15,19,23,25,5,11,13,17,21)]
data_h1_p2_ex <- data_h1_ex2[,c(2,4,8,10,16,20,24,26,6,12,14,18,22)]


#p_vals<-t(data_h1_p2_ex)
#heatmap_data <- t(data_h1_score2_ex)

p_vals<-data_h1_p2_ex
heatmap_data <- data_h1_score2_ex

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight =8,fontsize=8,border_color="black",scale="none",
                          fontsize_row=8, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77","#E77A77", "#ECBA84", "white","lightblue",'#336699',"blue")))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data >1.5, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))


##oligodendrocytes
oligo <- read.table("heatmap_oligodendrocytes.txt",header = TRUE,sep = "\t")


data_h1_ex <- oligo[,-1]
rownames(data_h1_ex) <- oligo[,1]


data_h1_ex2 <- filter(data_h1_ex,(ADHD_p < 0.05 & ADHD_score > 1.5)| (AN_p < 0.05 & AN_score > 1.5) |(AMD_p < 0.05 & AMD_score > 1.5)| (ASD_p < 0.05 & ASD_score > 1.5) | (DR_p < 0.05 & DR_score > 1.5) | (RD_p < 0.05 & RD_score > 1.5) | (BIP_p < 0.05 & BIP_score > 1.5) | (GLC_p < 0.05 & GLC_score  > 1.5) | (MDD_p < 0.05 & MDD_score > 1.5)| (SCZ_p < 0.05 & SCZ_score > 1.5) | (TS_p < 0.05 & TS_score > 1.5) | (OCD_p < 0.05 & OCD_score > 1.5)| (myopia_p < 0.05 & myopia_score > 1.5))

length(data_h1_ex2$ADHD_score)

#data_h1_score2_ex <- data_h1_ex2[,c(1,3,5,7,9,11,13,15,17,19,21,23,25)]
#data_h1_p2_ex <- data_h1_ex2[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
data_h1_score2_ex <- data_h1_ex2[,c(1,3,7,9,15,19,23,25,5,11,13,17,21)]
data_h1_p2_ex <- data_h1_ex2[,c(2,4,8,10,16,20,24,26,6,12,14,18,22)]


#p_vals<-t(data_h1_p2_ex)
#heatmap_data <- t(data_h1_score2_ex)

p_vals<-data_h1_p2_ex
heatmap_data <- data_h1_score2_ex

pearson_heatmap<-pheatmap(heatmap_data,cellwidth =8, cellheight =8,fontsize=8,border_color="black",scale="none",
                          fontsize_row=8, 
                          color=colorRampPalette(rev(c( "firebrick3","#E77A77", "#ECBA84", "white","lightblue",'#336699',"blue")))(102),
                          cluster_rows = T,
                          cluster_cols = F,
                          display_numbers = matrix(ifelse(p_vals <0.01 & heatmap_data >1.5, "**", ifelse(p_vals<0.05 & heatmap_data>1.5,"*","")), nrow(p_vals)))




####--------------Shared regulons and specific regulons across eye-brain diseases in a given cell types

##running for human brain single-cell multiomic data
#grn_outputs <- GRN_func(single_cell = pbmc_10x_brain,n_genes=10)

#all regulons
grn_outputs
length(grn_outputs$tf_names)
length(unique(grn_outputs$grn$regulons.meta.tf))
#extract regulons from grn_outputs for excitatory neurons

tf_left_shared <- c("BCL11A","STAT4","EGR4","GLIS1","SATB2","MEF2C","NPAS4",
                    "SMAD3","THRB","POU2F2","PKNOX2","TSHZ3","HLF","MYT1L","AEBP2")
tf_left_specific <- c("TBR1","ONECUT2","EMX1","ARNTL2","TRERF1","ZNF391","HES7",
                    "ZNF697","DACH2","TOX3","UBTF","NEUROD2","ZFPM2")

length(tf_left_shared)
length(tf_left_specific)
length(tf_left_shared2)
length(tf_left_specific2)

#extracting the gene list of each regulon for shared
M1 <- grn_outputs$grn[which(grn_outputs$grn$regulons.meta.tf == tf_left_shared[1]),]
M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
M1_regulon <- unique(M1_regulon)
M1_regulon
length(unique(M1_regulon))

write.csv(M1_regulon, file = "regulons_BCL11A_across13eye-brain_disease_excitatory.csv")


M1 <- grn_outputs$grn[which(grn_outputs$grn$regulons.meta.tf == tf_left_shared[2]),]
M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
M1_regulon <- unique(M1_regulon)
M1_regulon
length(unique(M1_regulon))

write.csv(M1_regulon, file = "regulons_STAT4_across13eye-brain_disease_excitatory.csv")



M1 <- grn_outputs$grn[which(grn_outputs$grn$regulons.meta.tf %in% tf_left_shared),]
M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
M1_regulon <- unique(M1_regulon)
M1_regulon
length(unique(M1_regulon))

write.csv(M1_regulon, file = "regulons_shared_across13eye-brain_disease_excitatory.csv")

M1 <- grn_outputs$grn[which(grn_outputs$grn$regulons.meta.tf %in% tf_left_specific),]
M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
M1_regulon <- unique(M1_regulon)
M1_regulon
length(unique(M1_regulon))
write.csv(M1_regulon, file = "regulons_specific_across13eye-brain_disease_excitatory.csv")

###


#extract regulons from grn_outputs for Microglia

tf_left_shared2 <- c("PRDM1","FOXP2","MEF2C","KLF6","SPI1","FLI1","JUNB",
                    "ARID3A","PLAG1","SMAD3","IRF8","NFKB2","ARID5A","TFEC","NFATC2","ETS2")
tf_left_specific2 <- c("BNC2","TAL1","JAZF1","RUNX2","KLF2","REL","HIVEP3",
                      "MAF","NR4A1","FOSB","SMAD7","ZNF217","CEBPA","FOSL1","ELF1",
                      "ZFHX3","NR4A2","RREB1")

M1 <- grn_outputs$grn[which(grn_outputs$grn$regulons.meta.tf == tf_left_shared2[2]),]
M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
M1_regulon <- unique(M1_regulon)
M1_regulon
length(unique(M1_regulon))

write.csv(M1_regulon, file = "regulons_FOXP2_across13eye-brain_disease_microglia.csv")


M1 <- grn_outputs$grn[which(grn_outputs$grn$regulons.meta.tf %in% tf_left_shared2),]
M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
M1_regulon <- unique(M1_regulon)
M1_regulon
length(unique(M1_regulon))
write.csv(M1_regulon, file = "regulons_shared_across13eye-brain_disease_microglia.csv")

M1 <- grn_outputs$grn[which(grn_outputs$grn$regulons.meta.tf %in% tf_left_specific2),]
M1_regulon <- c(unique(M1$regulons.meta.tf),M1$regulons.meta.target)
M1_regulon <- unique(M1_regulon)
M1_regulon
length(unique(M1_regulon))
write.csv(M1_regulon, file = "regulons_specific_across13eye-brain_disease_microliga.csv")

