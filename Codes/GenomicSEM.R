library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyverse)
library(ggtext)
library(wesanderson)
library(tidyr)
library(reshape2)  
#Downlaod and install GenomicSEM
install.packages("devtools")
library(devtools)
install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
setwd("/share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/")
datapath="/share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/"
gwas_list<-list.files(path=datapath,pattern = ".sumstats.gz")
#running LDSC
traits<-gwas_list
sample.prev<-c(.172,.234,.397,.101,.346,.276,.337,.407,.088,.027,.115,.436,.054)  #(样本患病率)
population.prev<-c(.05,.01,.01,.01,.15,.025,.008,.01,.123,.0176,.0251,.306,.0002) #(人口患病率)
ld<-"/share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/"
wld<-"/share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/"
trait.names<-c("ADHD","AN","ASD","BIP","MDD","OCD","TS","SCZ","AMD","DR","GLC","Myopia","RD")
LDSCoutput<-ldsc(traits=traits,sample.prev=sample.prev,population.prev=population.prev,ld=ld,wld=wld,trait.names=trait.names)
save(LDSCoutput,file="LDSCoutput_13.RData")
#LDSCoutput是包含5个命名变量的列表
#LDSCoutput$S是协方差矩阵。
#LDSCoutput$V这是 lavaan 预期格式的采样协方差矩阵。
#LDSCoutput$I是LDSC截距和交叉特征（即二元）截距的矩阵。
#LDSCoutput$N包含遗传力的样本数量 （N） 和共遗传力的 sqrt（N1N2）。
#LDSCoutput$m是用于构建 LD 分数的 SNP 数量。
library(Matrix)
library(stats)
load("LDSCoutput_13.RData")
#DWLS 最小二乘  ML 最大似然
#共因数分析
CommonFactor_DWLS <- commonfactor(covstruc = LDSCoutput, estimation="DWLS")
CommonFactor_DWLS
CommonFactor_ML <- commonfactor(covstruc = LDSCoutput, estimation="ML")
CommonFactor_ML
write.table(CommonFactor_DWLS$results, file ="CommonFactor_DWLS.txt", sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(CommonFactor_ML$results, file = "CommonFactor_ML.txt", sep = "\t",col.names = TRUE,row.names = FALSE)
#探索性因素分析
#smooth the S matrix for EFA using the nearPD function in the Matrix package. 
require(Matrix)
Ssmooth <- as.matrix((nearPD(LDSCoutput$S, corr = FALSE))$mat)
Ssmooth

#碎石图
install.packages("nFactors")
library(nFactors)
ev <- eigen(cor(Ssmooth)) # 获取特征值
ap <- parallel(subject=nrow(Ssmooth),var=ncol(Ssmooth),
               rep=100,cent=.05) # subject指样本个数，var是指变量个数
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)# 确定探索性因子分析中应保留的因子
plotnScree(nS) # 绘制碎石图
#根据碎石图以及相关知识 将因素分为五组
#路径图
require(stats)
EFA <- factanal(covmat = Ssmooth, factors = 5, rotation = "promax")
EFA$loadings
fa.diagram(EFA$loadings,simple = F)



#验证性分析

CFAofEFA <- 'F1 =~ AMD + DR + GLC +TS
             F2 =~ ADHD + ASD + MDD
             F3 =~ OCD + TS + AN + GLC
             F4 =~ BIP + MDD + SCZ
             F5 =~ Myopia + RD + ASD
F1~~F2
F1~~F3
F1~~F4
F1~~F5
F2~~F5
F2~~F4
F2~~F3
F3~~F4
F3~~F5
F4~~F5
F1 ~~ 1*F1
F2 ~~ 1*F2
F3 ~~ 1*F3
Myopia ~~ a*Myopia
a > .001
'


#run the model
Anthro_DWLS <- usermodel(LDSCoutput, estimation = "DWLS", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
Anthro_ML <- usermodel(LDSCoutput, estimation = "ML", model = CFAofEFA, CFIcalc = TRUE, std.lv = TRUE, imp_cov = FALSE)
#print the results
Anthro_DWLS
Anthro_ML
write.table(Anthro_DWLS$results, file = "CFA_DWLS.txt", sep = "\t",col.names = TRUE,row.names = FALSE)
write.table(Anthro_ML$results, file = "CFA_DWLS.txt", sep = "\t",col.names = TRUE,row.names = FALSE)
