#!/usr/bin/env bash

#MAGMA_gene-based association analysis for 13 eye-brain diseases


#DIRECTORY
export MAGMA_DIR=/share/pub/mayl/Sherlock/MAGMA
export DATA=/share/pub/mayl/Sherlock/MAGMA_test
export OUTPUT=/share/pub/mayl/Sherlock/MAGMA_test

#Formating
#cd $DATA

#generating a location file including three Columns: SNP, CHR, POS
#  nohup gawk '{print $1, $7, $8}'  ../GWAS_UKBiobank_summary_final   > GWAS_UKBiobank_summary_final.hg19.location &

#generating a --pval file including two Columns: SNP, P
#If MAGMA detects a header in the file it will look for SNP IDs and p-values in the SNP and P column respectively. 
#If no header is found it will use the first column for SNP IDs and the second column for p-values.
#   nohup gawk '{print $1, $4}' ../GWAS_UKBiobank_summary_final > GWAS_UKBiobank_summary_final.results_Pval &


#MAGMA annotation:

$MAGMA_DIR/magma \
    --snp-loc  $DATA/GWAS_UKBiobank_summary_final.hg19.location  \
    --annotate window=20,20 --gene-loc $MAGMA_DIR/NCBI37.3.gene.loc \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation  


#gene-based association analysi:
$MAGMA_DIR/magma \
    --bfile $MAGMA_DIR/1000G_data/g1000_eur \
    --pval $DATA/GWAS_UKBiobank_summary_final.results_Pval \
    N=13239 \
    --gene-annot   $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_annotation.genes.annot  \
    --out $OUTPUT/GWAS_UKBiobank_summary_final.hg19_SNP_Gene_Analysis_P
