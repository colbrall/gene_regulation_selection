#!/bin/bash

#make sure R >3.3 is available
# source activate R_peer

chr=$1
study=$2
expr_RDS=$3
geno=$4
gene_annot=$5
snp_annot=$6
n_k_folds=$7
alpha=$8
out_dir=$9
snpset=${10}
window=${11}
home_dir=${12}

gen=$geno$chr'.txt'
snps=$snp_annot$chr'.RDS'

Rscript ./create_model.R $study $expr_RDS $gen $gene_annot $snps \
    $n_k_folds $alpha $out_dir $chr $snpset $window
