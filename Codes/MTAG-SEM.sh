#mtag下载
conda create -n py27 python=2.7 #创建py27
conda activate py27 #激活
conda install numpy
conda install scipy 
conda install pandas
conda install argparse
conda install bitarray
conda install joblib
conda install libgfortran==1
cd /share/pub1/gouxx/gouxx/softwares/
git clone https://github.com/omeed-maghzian/mtag.git
python mtag.py -h
#安装完成
#MTAG必须要Z值输入，beta/or输入已在后续版本中移除
#输入文件包含以下几列：snpid, chr, bpos, a1, a2, freq, z, pval 和 n；
#snpid指SNP的ID,一般用RS表示；
#chr指染色体；
#bpos指SNP的位置；
#a1指效应位点；
#a2指非效应位点；
#freq指a1的频率;
#z指zscore，可通过beta/se获得;
#pval指p值；
#n指有效样本数；


#ADHD ASD BIP MDD OCD TS 数据格式为 or 无zscore
#AN SCZ 数据格式为 beta 无zscore
#AMD DR GLC RD Myopia 数据格式为 beta

#beta=log（or）  zscore=beta/se  zscore=(log(or))/se

#Rstudio中处理
#计算补齐缺失的zscore  z=log(or)/se   z=beta/se




#MTAG运行 按GemomicSEM分组结果运行mtag
cd /share/pub1/gouxx/gouxx/GWAS/brain_eye/
/share/pub1/gouxx/gouxx/softwares/mtag/mtag.py 
--sumstats ldsc_eye_AMD.vcf,ldsc_eye_RD.vcf,ldsc_eye_GLC.vcf,Z_ldsc_brain_TS.vcf \
--out /share/pub1/gouxx/gouxx/temp/GenomicSEMgroups/factor1 \
--n_min 0.0 \
--stream_stdout &

/share/pub1/gouxx/gouxx/softwares/mtag/mtag.py 
--sumstats Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_MDD.vcf \
--out /share/pub1/gouxx/gouxx/temp/GenomicSEMgroups/factor2 \
--n_min 0.0 \
--stream_stdout &

/share/pub1/gouxx/gouxx/softwares/mtag/mtag.py 
--sumstats ldsc_eye_AMD.vcf,ldsc_eye_RD.vcf,ldsc_eye_GLC.vcf,Z_ldsc_brain_TS.vcf \
--out /share/pub1/gouxx/gouxx/temp/GenomicSEMgroups/factor3 \
--n_min 0.0 \
--stream_stdout &

/share/pub1/gouxx/gouxx/softwares/mtag/mtag.py 
--sumstats Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_SCZ.vcf \
--out /share/pub1/gouxx/gouxx/temp/GenomicSEMgroups/factor4 \
--n_min 0.0 \
--stream_stdout &

/share/pub1/gouxx/gouxx/softwares/mtag/mtag.py 
--sumstats Z_ldsc_brain_ASD.vcf,ldsc_eye_myopia.vcf,ldsc_eye_RD.vcf \
--out /share/pub1/gouxx/gouxx/temp/GenomicSEMgroups/factor5 \
--n_min 0.0 \
--stream_stdout &


