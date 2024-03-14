#!/bin/bash
#PBS -N ldsc-rg
#PBS -q workq
#PBS -j oe
#PBS -l ncpus=1
#PBS -l nodes=node03

source /share/pub1/gouxx/gouxx/miniconda3/bin/activate ldsc
cd /share/pub1/gouxx/gouxx/GWAS/brain_eye_sumstats/

for fileeye in ldsc_Eye_*.sumstats.gz
do
	filename=$(basename "$fileeye" .sumstats.gz)
  for filebrain in ldsc_*.sumstats.gz
  do
    fileeye="$fileeye,$filebrain"
  done
  python /share/pub1/gouxx/gouxx/softwares/ldsc/ldsc.py \
  --rg $fileeye \
  --ref-ld-chr /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/ \
  --w-ld-chr /share2/pub/jiangdp/jiangdp/ldsc-master/eur_w_ld_chr/ \
  --out /share/pub1/gouxx/gouxx/outcomes/ldsc-brain-eye/$filename
done
