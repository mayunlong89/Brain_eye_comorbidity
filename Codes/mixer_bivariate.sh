#!/usr/bin/env bash
#PBS -N x_bi_ADHD_AMD19_20.pbs
#PBS -q workq
#PBS -l ncpus=4
#PBS -l mem=50gb  

export TOOLS=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/mixer
export python_convert=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/LDSC_file/python_convert
export OUTPUTS=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/LDSC_file/mixer_fit_models
export RESOURCES=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/LDSC_file/1000G_EUR_Phase3_plink

##This is the vital version R for using
conda activate /share/pub/mayl/00_conda_source/R4.2.3
#source activate /share/pub/mayl/00_conda_source/R4.2.3


trait1=ADHD
trait2=AMD

for i in $(seq 19 20)  
do   
SLURM_ARRAY_TASK_ID=$i
#@@2 Bivariate (cross-trait) analysis 
#Fit the model: 
python3 $TOOLS/precimed/mixer.py fit2 --trait1-file $python_convert/PGC_${trait1}.EUR_qc_noMHC.csv.gz --trait2-file $python_convert/PGC_${trait2}.EUR_qc_noMHC.csv.gz --trait1-params-file $OUTPUTS/PGC_${trait1}.EUR_qc_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json --trait2-params-file $OUTPUTS/PGC_${trait2}.EUR_qc_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json --out $OUTPUTS/Bivariate/PGC_${trait1}.EUR_qc_noMHC_vs_PGC_${trait2}.EUR_qc_noMHC.fit.rep${SLURM_ARRAY_TASK_ID} --extract $RESOURCES/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps --bim-file $RESOURCES/1000G.EUR.QC.@.bim --ld-file $RESOURCES/1000G.EUR.QC.@.run4.ld --lib  $TOOLS/src/build/lib/libbgmg.so  

#Apply the model to the entire set of SNPs, without constraining to LDSR/w_hm3.justrs: 
python3 $TOOLS/precimed/mixer.py test2 --trait1-file $python_convert/PGC_${trait1}.EUR_qc_noMHC.csv.gz --trait2-file $python_convert/PGC_${trait2}.EUR_qc_noMHC.csv.gz --load-params-file $OUTPUTS/Bivariate/PGC_${trait1}.EUR_qc_noMHC_vs_PGC_${trait2}.EUR_qc_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json --out $OUTPUTS/Bivariate/PGC_${trait1}.EUR_qc_noMHC_vs_PGC_${trait2}.EUR_qc_noMHC.test.rep${SLURM_ARRAY_TASK_ID} --bim-file $RESOURCES/1000G.EUR.QC.@.bim --ld-file $RESOURCES/1000G.EUR.QC.@.run4.ld --lib  $TOOLS/src/build/lib/libbgmg.so


done
