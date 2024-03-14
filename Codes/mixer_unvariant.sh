#!/usr/bin/env bash
#PBS -N x_mix_ASD.pbs
#PBS -q workq
#PBS -l ncpus=4
#PBS -l mem=40gb  


export TOOLS=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/mixer
export python_convert=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/LDSC_file/python_convert
#export OUTPUTS=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/LDSC_file/mixer_fit_models
export OUTPUTS=/share2/pub/zhenggw/zhenggw/mixer_test
export RESOURCES=/share/pub/mayl/eye-brain/MiXeR/boost_1_84_0/LDSC_file/1000G_EUR_Phase3_plink

##This is the vital version R for using
#conda activate /share/pub/mayl/00_conda_source/R4.2.3
source activate /share/pub/mayl/00_conda_source/R4.2.3

name=ASD

for i in $(seq 1 19)  
do   
SLURM_ARRAY_TASK_ID=$i

##-Univariate analysis
#Fit the model:
python3 $TOOLS/precimed/mixer.py fit1 --trait1-file $python_convert/PGC_${name}.EUR_qc_noMHC.csv.gz --out $OUTPUTS/PGC_${name}.EUR_qc_noMHC.fit.rep${SLURM_ARRAY_TASK_ID} --extract $RESOURCES/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${SLURM_ARRAY_TASK_ID}.snps --bim-file $RESOURCES/1000G.EUR.QC.@.bim --ld-file $RESOURCES/1000G.EUR.QC.@.run4.ld --lib  $TOOLS/src/build/lib/libbgmg.so  

#Apply the model to the entire set of SNPs, without constraining to LDSR/w_hm3.justrs:
python3 $TOOLS/precimed/mixer.py test1 --trait1-file $python_convert/PGC_${name}.EUR_qc_noMHC.csv.gz --load-params-file $OUTPUTS/PGC_${name}.EUR_qc_noMHC.fit.rep${SLURM_ARRAY_TASK_ID}.json --out $OUTPUTS/PGC_${name}.EUR_qc_noMHC.test.rep${SLURM_ARRAY_TASK_ID} --bim-file $RESOURCES/1000G.EUR.QC.@.bim --ld-file $RESOURCES/1000G.EUR.QC.@.run4.ld --lib $TOOLS/src/build/lib/libbgmg.so 

done
