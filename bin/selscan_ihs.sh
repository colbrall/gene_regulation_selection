#!/bin/bash
# adapted from Lin's script, calls her ancestral_match.R

pop=$1

data=/project/mathilab/data/1000g/NYGC_phased
annots=/project/mathilab/data/1000g/NYGC
chain=/project/mathilab/ppoyraz/lift/hg38ToPanTro6.over.chain.gz

less ../data/1kG_pops/1kG_populations.txt | grep ${pop} |cut -f 1 > tmp/${pop}_subset.txt

for chr in {1..22}; do
# separate target population
  bcftools view ${data}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz -S tmp/${pop}_subset.txt -m2 -M2 -v snps |bcftools view - -g ^miss -q 0.05 | bgzip > tmp/iHS/1kG_chr${chr}_${pop}.vcf_not_annot.gz;
  bcftools index tmp/iHS/1kG_chr${chr}_${pop}.vcf_not_annot.gz
  bcftools annotate -Oz -a ${annots}/CCDG_13607_B01_GRM_WGS_2019-02-19_chr${chr}.recalibrated_variants.vcf.gz -c ID tmp/iHS/1kG_chr${chr}_${pop}.vcf_not_annot.gz  > tmp/iHS/1kG_chr${chr}_${pop}.vcf.gz
  # rm tmp/iHS/1kG_chr${chr}_${pop}.vcf_not_annot.gz*;


  zcat tmp/iHS/1kG_chr${chr}_${pop}.vcf.gz | awk -v OFS='\t' -v FS='\t' '/^[^#]/ {print $1, $2, $2, $3}' > tmp/liftover/chr${chr}_${pop}.bed
  ./../bin/liftOver tmp/liftover/chr${chr}_${pop}.bed ${chain} tmp/liftover/panTro_chr${chr}_${pop}.bed tmp/liftover/unlifted_chr${chr}_${pop}.bed

# polarize
  Rscript ./bin/ancestral_match.R tmp/liftover/panTro_chr${chr}_${pop}.bed tmp/liftover/unlifted_chr${chr}_${pop}.bed tmp/iHS/1kG_chr${chr}_${pop}.vcf.gz ${pop}
  awk '{print $1}' tmp/chr${chr}_${pop}_AA.txt | grep ^rs > tmp/liftover/chr${chr}_${pop}_AA_filter.txt
  plink2 --vcf tmp/iHS/1kG_chr${chr}_${pop}.vcf.gz --extract tmp/liftover/chr${chr}_${pop}_AA_filter.txt --ref-allele force tmp/chr${chr}_${pop}_AA.txt --recode vcf --out tmp/iHS/chr${chr}_${pop}_AA_filtered
  bgzip -f tmp/iHS/chr${chr}_${pop}_AA_filtered.vcf
  # rm tmp/iHS/1kG_chr${chr}_${pop}.vcf.gz

# interpolate rr map, run selscan
  zcat tmp/iHS/chr${chr}_${pop}_AA_filtered.vcf.gz | grep -v "#" |cut -f 2 > tmp/${pop}_chr${chr}_newpos.txt;
  ./bin/interpolate_map.R chr${chr} ../../data/maps/1kG_hg38/${pop}/chr${chr}_hg38.map tmp/${pop}_chr${chr}_newpos.txt tmp/interp_maps/${pop}_chr${chr}.map
  ./bin/selscan/selscan --ihs --map tmp/interp_maps/${pop}_chr${chr}.map --vcf tmp/iHS/chr${chr}_${pop}_AA_filtered.vcf.gz --out results/selection/iHS/${pop}_chr${chr}
  # rm tmp/${pop}_chr${chr}_newpos.txt
  # rm tmp/iHS/chr${chr}_${pop}_AA_filtered.vcf.gz
done
