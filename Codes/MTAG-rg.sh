#!/bin/bash
#PBS -N mtag
#PBS -q workq
#PBS -j oe
#PBS -l ncpus=10
#PBS -l mem=100g


source /share/pub1/gouxx/gouxx/miniconda3/bin/activate py27
cd /share/pub1/gouxx/gouxx/GWAS/mtag_brain_eye/
python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_ASD.vcf,Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_SCZ.vcf,Z_ldsc_brain_TS.vcf,ldsc_eye_AMD.vcf,ldsc_eye_DR.vcf,ldsc_eye_myopia.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/ADHD \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_AN.vcf,Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_OCD.vcf,Z_ldsc_brain_SCZ.vcf,ldsc_eye_DR.vcf,ldsc_eye_myopia.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/AN \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_ASD.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_SCZ.vcf,Z_ldsc_brain_TS.vcf,ldsc_eye_myopia.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/ASD \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_ASD.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_OCD.vcf,Z_ldsc_brain_SCZ.vcf,ldsc_eye_AMD.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/BIP \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_OCD.vcf,Z_ldsc_brain_SCZ.vcf,Z_ldsc_brain_TS.vcf,ldsc_eye_AMD.vcf,ldsc_eye_DR.vcf,ldsc_eye_GLC.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/MDD \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_OCD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_SCZ.vcf,Z_ldsc_brain_TS.vcf,ldsc_eye_DR.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/OCD \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_SCZ.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_ASD.vcf,Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_OCD.vcf,Z_ldsc_brain_TS.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/SCZ \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats Z_ldsc_brain_TS.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_ASD.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_OCD.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/TS \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats ldsc_eye_AMD.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_BIP.vcf,Z_ldsc_brain_MDD.vcf,ldsc_eye_GLC.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/AMD \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats ldsc_eye_DR.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_MDD.vcf,Z_ldsc_brain_OCD.vcf,ldsc_eye_GLC.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/DR \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats ldsc_eye_GLC.vcf,Z_ldsc_brain_MDD.vcf,ldsc_eye_AMD.vcf,ldsc_eye_DR.vcf,ldsc_eye_RD.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/GLC \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats ldsc_eye_myopia.vcf,Z_ldsc_brain_ADHD.vcf,Z_ldsc_brain_AN.vcf,Z_ldsc_brain_ASD.vcf,ldsc_eye_RD.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/myopia \
--n_min 0.0 \
--stream_stdout 

python /share/pub1/gouxx/gouxx/softwares/mtag/mtag.py \
--sumstats ldsc_eye_RD.vcf,ldsc_eye_GLC.vcf,ldsc_eye_myopia.vcf \
--out /share/pub1/gouxx/gouxx/outcomes/RD \
--n_min 0.0 \
--stream_stdout 
