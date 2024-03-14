setwd("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\")


##Genetic Correlations between LDSC and GENOVA
data_cor_sig <- read.table("./01-GNOVA/LDSC_GNOVA_cor_sig.txt",header = TRUE)
data_cor_all <- read.table("./01-GNOVA/LDSC_GNOVA_cor_all.txt",header = TRUE)


cor.test(data_cor_sig$LDSC,data_cor_sig$GNOVA)
cor.test(data_cor_all$LSDC,data_cor_all$GNOVA)

plot(data_cor_sig$LDSC~data_cor_sig$GNOVA,col="red")
abline(lm(data_cor_sig$LDSC~data_cor_sig$GNOVA))

plot(data_cor_all$LSDC~data_cor_all$GNOVA)

 
library(ggplot2) ##s-------------------------------------------Scatter plot-smooth line
#library(ggpmisc) #install.packages("ggpmisc")
ggplot(data=data_cor_sig,mapping = aes(x=LDSC,y=GNOVA))+
  geom_point(pch=19,col= "#69B9DA",alpha=1,cex=3.0)+  #alpha 设置透明度，越小越透明
  geom_smooth(method = "lm",formula = y~x, col="red",alpha=0.2)+
  theme_classic()  ###这个放到前面不要放到最后，否则转角angle不起作用


library(ggplot2) ##s-------------------------------------------Scatter plot-smooth line
#library(ggpmisc) #install.packages("ggpmisc")
ggplot(data=data_cor_all,mapping = aes(x=LSDC,y=GNOVA))+
  geom_point(pch=19,col= "#69B9DA",alpha=1,cex=3.0)+  #alpha 设置透明度，越小越透明
  geom_smooth(method = "lm",formula = y~x, col="red",alpha=0.2)+
  theme_classic()  ###这个放到前面不要放到最后，否则转角angle不起作用
