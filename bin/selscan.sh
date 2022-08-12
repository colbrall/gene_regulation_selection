#!/bin/bash

pop=$1
less ../data/1kG_pops/1kG_populations.txt | grep ${pop} |cut -f 1 > tmp/${pop}_subset.txt
echo ${pop}

for chr in {21..22}; do
echo ${chr};
bcftools view ../../data/1000g/NYGC/CCDG_13607_B01_GRM_WGS_2019-02-19_chr${chr}.recalibrated_variants.vcf.gz -S tmp/${pop}_subset.txt -m2 -M2 -v snps | bcftools view - -g ^miss -q 0.05 | bgzip > tmp/1kG_chr${chr}_${pop}.vcf.gz;
./bin/selscan/selscan --nsl --vcf tmp/1kG_chr${chr}_${pop}.vcf.gz --out results/selection/nSL/${pop}_chr${chr} --unphased;
rm tmp/1kG_chr${chr}_${pop}.vcf.gz;
done
