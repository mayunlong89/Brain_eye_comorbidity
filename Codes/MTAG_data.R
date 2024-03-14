#计算补齐缺失的zscore  z=-log(or)/se   z=beta/se
setwd('/share/pub1/gouxx/gouxx/GWAS/mtag_brain_eye')
#ADHD<-read.table('ldsc_brain_ADHD.vcf',header = T)
#ADHD
#a<-log(ADHD$or)
#a
#z<-a/ADHD$se
#b
#Z_ADHD<-cbind(ADHD,z)



filenames<-list.files(pattern=".vcf")
file.list<-list()
length(file.list)<-length(filenames)
names(file.list)<-filenames
for (i in 1:length(filenames)){
  file.list[[i]]<-read.table(filenames[i],header = T)
}
#AN SCZ 为beta  file.list[2&7]  计算z
file.list[[1]]$z<-log(file.list[[1]]$or)/file.list[[1]]$se
file.list[[2]]$z<-file.list[[2]]$beta/file.list[[2]]$se
file.list[[3]]$z<-log(file.list[[3]]$or)/file.list[[3]]$se
file.list[[4]]$z<-log(file.list[[4]]$or)/file.list[[4]]$se
file.list[[5]]$z<-log(file.list[[5]]$or)/file.list[[5]]$se
file.list[[6]]$z<-log(file.list[[6]]$or)/file.list[[6]]$se
file.list[[7]]$z<-file.list[[7]]$beta/file.list[[7]]$se
file.list[[8]]$z<-log(file.list[[8]]$or)/file.list[[8]]$se

#将列表转化为数据框
#x<-data.frame(file.list[1])
#colnames(x)<-colnames(file.list[[1]])
Z_ldsc_brain_ADHD<-data.frame(file.list[1])
colnames(Z_ldsc_brain_ADHD)<-colnames(file.list[[1]])
Z_ldsc_brain_AN<-data.frame(file.list[2])
colnames(Z_ldsc_brain_AN)<-colnames(file.list[[2]])
Z_ldsc_brain_ASD<-data.frame(file.list[3])
colnames(Z_ldsc_brain_ASD)<-colnames(file.list[[3]])
Z_ldsc_brain_BIP<-data.frame(file.list[4])
colnames(Z_ldsc_brain_BIP)<-colnames(file.list[[4]])
Z_ldsc_brain_MDD<-data.frame(file.list[5])
colnames(Z_ldsc_brain_MDD)<-colnames(file.list[[5]])
Z_ldsc_brain_OCD<-data.frame(file.list[6])
colnames(Z_ldsc_brain_OCD)<-colnames(file.list[[6]])
Z_ldsc_brain_SCZ<-data.frame(file.list[7])
colnames(Z_ldsc_brain_SCZ)<-colnames(file.list[[7]])
Z_ldsc_brain_TS<-data.frame(file.list[8])
colnames(Z_ldsc_brain_TS)<-colnames(file.list[[8]])

write.table(Z_ldsc_brain_ADHD,file = "Z_ldsc_brain_ADHD.vcf",col.names = T,row.names = F,quote = F)
write.table(Z_ldsc_brain_AN,file = "Z_ldsc_brain_AN.vcf",col.names = T,row.names = F,quote = F)
write.table(Z_ldsc_brain_ASD,file = "Z_ldsc_brain_ASD.vcf",col.names = T,row.names = F,quote = F)
write.table(Z_ldsc_brain_BIP,file = "Z_ldsc_brain_BIP.vcf",col.names = T,row.names = F,quote = F)
write.table(Z_ldsc_brain_MDD,file = "Z_ldsc_brain_MDD.vcf",col.names = T,row.names = F,quote = F)
write.table(Z_ldsc_brain_OCD,file = "Z_ldsc_brain_OCD.vcf",col.names = T,row.names = F,quote = F)
write.table(Z_ldsc_brain_SCZ,file = "Z_ldsc_brain_SCZ.vcf",col.names = T,row.names = F,quote = F)
write.table(Z_ldsc_brain_TS,file = "Z_ldsc_brain_TS.vcf",col.names = T,row.names = F,quote = F)
