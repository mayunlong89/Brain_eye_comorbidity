##Sample sizes bias estimation for LDSC and GENOVA
data_sample_bias <- read.table("./01-GNOVA/sample_bias_test.txt",header = TRUE)
 

cor.test(data_sample_bias$sample_sizes,data_sample_bias$LDSC)
cor.test(data_sample_bias$sample_sizes,data_sample_bias$GENOVA)

library(ggplot2) ##s-------------------------------------------Scatter plot-smooth line
#library(ggpmisc) #install.packages("ggpmisc")
ggplot(data=data_sample_bias,mapping = aes(x=sample_sizes,y=LDSC))+
  geom_point(pch=19,col= "#69B9DA",alpha=1,cex=3.0)+  #alpha 设置透明度，越小越透明
  #geom_smooth(method = "lm",formula = y~x, col="red",alpha=0.2)+
  geom_text(aes(sample_sizes, LDSC, label=Trait_pairs))+
  theme_classic()  ###这个放到前面不要放到最后，否则转角angle不起作用


library(ggplot2) ##s-------------------------------------------Scatter plot-smooth line
#library(ggpmisc) #install.packages("ggpmisc")
ggplot(data=data_sample_bias,mapping = aes(x=sample_sizes,y=GENOVA))+
  geom_point(pch=19,col= "#69B9DA",alpha=1,cex=3.0)+  #alpha 设置透明度，越小越透明
  #geom_smooth(method = "lm",formula = y~x, col="red",alpha=0.2)+
  geom_text(aes(sample_sizes, GENOVA, label=Trait_pairs))+
  theme_classic()  ###这个放到前面不要放到最后，否则转角angle不起作用
