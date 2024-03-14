#ldsc数据格式
#数据地址 cd /share2/pub/jiangdp/jiangdp/GWAS/Xuezb_EYE_Cell_Brain/raw_data/
cd /share/pub1/gouxx/gouxx/GWAS/brain_eye_raw/
#ldsc sumstats格式 需要下载ldsc
conda activate ldsc

munge_sumstats.py \
--sumstats ldsc_brain_SCZ.vcf \
--chunksize 500000 \
--N 127906 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_SCZ

munge_sumstats.py \
--sumstats ldsc_brain_ADHD.vcf \
--N 225534 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_ADHD \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_brain_BIP.vcf \
--N 413466 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_BIP \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_brain_MDD.vcf \
--N 173005 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_MDD \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_brain_OCD.vcf \
--N 9725 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_OCD \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_brain_TS.vcf \
--N 14307 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_TS \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_brain_AN.vcf \
--N 72517 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_AN \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_brain_ASD.vcf \
--N 46351 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Brain_ASD \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_eye_AMD.vcf \
--ignore z \
--N 66387 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Eye_AMD \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_eye_DR.vcf \
--ignore z \
--N 62229 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Eye_DR \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_eye_GLC.vcf \
--ignore z \
--N 68390 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Eye_GLC \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_eye_myopia.vcf \
--ignore z \
--N 64268 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Eye_Myopia \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist

munge_sumstats.py \
--sumstats ldsc_eye_RD.vcf \
--ignore z \
--N 64011 \
--out /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/ldsc_Eye_RD \
--chunksize 500000 \
--merge-alleles /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/w_hm3.snplist
