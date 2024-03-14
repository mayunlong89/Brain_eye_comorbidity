#读取sum_brain_eye.txt 
setwd("/share/pub1/gouxx/gouxx/outcomes/ldsc-brain-eye/")
ldsc_sum1<-read.table("sum_brain_eye.txt",header = T)
p1<-rep(c('ADHD','AN','ASD','BIP','MDD','OCD','SCZ','TS','AMD','DR','GCL','Myopia','RD'),each=13)
p2<-rep(c('ADHD','AN','ASD','BIP','MDD','OCD','SCZ','TS','AMD','DR','GCL','Myopia','RD'),times=13)
ldsc_sum1$p1<-p1
ldsc_sum1$p2<-p2
#绘制相关性图
library(corrplot)
corr_data<-ldsc_sum1[,c('p1','p2','rg')]
library(tidyr)
corr_wide<-spread(corr_data,p2,rg)
row.names(corr_wide)<-corr_wide$p1
library(dplyr)
corr_wide<-corr_wide %>% select(-p1)
ldsc_rg<-corrplot(as.matrix(corr_wide), is.corr = FALSE, method = 'circle', type = 'upper',
         tl.col="black")
ggsave(ldsc_rg,filename = "ldsc_rg.pdf")
#绘制桑葚图
#按p<0.05过滤
setwd("/share/pub1/gouxx/gouxx/outcomes/ldsc-brain-eye/")
ldsc_sum2<-read.table("sum_brain.txt",header = T)
p11<-rep(c('ADHD','AN','ASD','BIP','MDD','OCD','SCZ','TS'),each=13)
p22<-rep(c('ADHD','AN','ASD','BIP','MDD','OCD','SCZ','TS','AMD','DR','GCL','Myopia','RD'),times=8)
ldsc_sum2$p1<-p11
ldsc_sum2$p2<-p22
sk_data <- as.data.frame(ldsc_sum2[ldsc_sum2$p<0.05,c("p1","p2","rg")])
sk_data$rg <- abs(sk_data$rg)
sk_data<-subset(sk_data,sk_data$p2=='AMD' | sk_data$p2=='GCL' | sk_data$p2=='DR' | sk_data$p2=='Myopia')

library(ggplot2)
library(ggalluvial)

is_alluvia_form(as.data.frame(sk_data), axes = 1:2, silent = TRUE)
ldsc_sankey<-ggplot(as.data.frame(sk_data),
       aes(y = rg, axis1 =p1, axis2 = p2)) +
  geom_alluvium(width = 1/6, aes(fill = p1),curve_type = "sigmoid") +
  geom_stratum(width = 1/6, color='white', aes(fill=p2)) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_fill_brewer(type="qual",palette = "Set3",guide=FALSE) +
  theme_void()
ggsave(ldsc_sankey,filename = "ldsc_sankey.pdf")
