#!/usr/bin/env python
# @author Laura Colbran, 10/2020
# adapted from train_models.py to run on bsub with Python 3

import subprocess
import time

from model_parameters_1240k_Zhang import *

CMD = 'bsub -W 24:00 -M 22000 -e models_{0}_%J.error -o models_{0}_%J.out -J models_{0} ' + \
    '"./train_predixcan.sh {0} {1} {2} {3} {4} {5} ' + \
    '{6} {7} {8} {9} {10} {11}"'

geno_prefixes = list(GENOTYPE_INTER_PREFIX)
gene_annot = INTER_DIR + GENE_ANNOT_INTER2
for i, study in enumerate(STUDY_NAMES):
    # for j in range(1,23):
    for j in range(15,16):
        expression_RDS = INTER_DIR + EXPRESSION_INTER_DIR + study + ".RDS"
        geno = INTER_DIR + GENOTYPE_INTER_DIR + geno_prefixes[0] + '.chr'
        snp_annot = INTER_DIR + SNP_ANN_INTER_DIR + SNP_ANN_INTER_PREFIX2
        cmd = CMD.format(j,study,expression_RDS,geno,gene_annot,snp_annot,
            N_K_FOLDS,ALPHA,MODEL_BY_CHR_DIR,SNPSET,WINDOW,HOME_DIR)
        print(cmd)
        subprocess.call(cmd, shell=True)
        time.sleep(2)
