#!/usr/bin/env python
# updated to run on Slurm with Python 3 by Laura Colbran 12/2018
# ** make sure the version of R you're running is earlier than 3.6-- struggles with readRDS later in pipeline if it's 3.6

import subprocess
import sys

from model_parameters_v8_1kGvariants import *

# Process gene annotation----------------------------------------------
print("Parsing gene annotation...")
subprocess.call(
    ['../data_processing/parse_gtf_v8.py',
    GENE_ANN_DIR + GENE_ANNOTATION_FN,
    INTER_DIR + GENE_ANNOT_INTER1
    ])
print("Saving each gene annotation file as RDS object")
subprocess.call(
    ['Rscript', '../data_processing/geno_annot_to_RDS.R',
    INTER_DIR + GENE_ANNOT_INTER1,
    INTER_DIR + GENE_ANNOT_INTER2])

# Process snp annotation-----------------------------------------------
print("Splitting SNP annotation file up by chromosome...")
subprocess.call(
  ['../data_processing/split_snp_annot_by_chr_v8.py',
  SNP_ANN_DIR + SNP_ANNOTATION_FN,
  INTER_DIR + SNP_ANN_INTER_DIR + SNP_ANN_INTER_PREFIX1
  ])
print("Saving each snp annotation file as RDS object")
subprocess.call(
  ['Rscript', '../data_processing/snp_annot_to_RDS.R',
  INTER_DIR + SNP_ANN_INTER_DIR + SNP_ANN_INTER_PREFIX2])

# Process genotype files-----------------------------------------------/
print("Splitting genotype files up by chromosome...")
for i in range(len(GENOTYPE_FNS)):
    subprocess.call(
        ['../data_processing/split_genotype_by_chr_v8.py',
        GENOTYPE_INPUT_DIR + GENOTYPE_FNS[i],
        INTER_DIR + GENOTYPE_INTER_DIR + GENOTYPE_INTER_PREFIX[i]])

# Process expression files---------------------------------------------
print("Transposing expression data and saving as RDS object...")
for i in range(len(STUDY_NAMES)):
  subprocess.call(
      ['Rscript', '../data_processing/expr_to_transposed_RDS_v8.R',
      STUDY_NAMES[i],
      EXPRESSION_INPUT_DIR + STUDY_NAMES[i],
      INTER_DIR + EXPRESSION_INTER_DIR])

# Create metadata files------------------------------------------------/
print("Creating metadata files...")
geno_prefix = list(GENOTYPE_INTER_PREFIX)
for i in range(len(STUDY_NAMES)):
    command = ' '.join(['../data_processing/create_meta_data.py',
        '--geno', INTER_DIR + GENOTYPE_INTER_DIR + geno_prefix[0] + '.chr22.txt',
        '--expr', INTER_DIR + EXPRESSION_INTER_DIR + STUDY_NAMES[i] + ".txt",
        '--snpset', SNPSET,
        '--alpha', ALPHA,
        '--n_k_folds', N_K_FOLDS,
        '--rsid_label', RSID_LABEL,
        '--window', WINDOW,
        '--out_prefix', OUTPUT_DIR + 'allMetaData/' + STUDY_NAMES[i]])
    subprocess.call(command, shell=True)

print("Done!")
