#!/usr/bin/env python
# @author Laura Colbran, 12/2018
# adapted from train_models.py to run on Slurm with Python 3

import subprocess
import time

# from model_parameters import *
# from model_parameters_normalized import *
# from model_parameters_1240k import *
# from model_parameters_top500k import *
# from model_parameters_v8_all import *
from model_parameters_v8_1kGvariants import *

CMD = 'sbatch --export=study={0},expr_RDS={1},geno={2},gene_annot={3},snp_annot={4},' + \
    'n_k_folds={5},alpha={6},out_dir={7},snpset={8},window={9} ' + \
    '-J {0}_models -D {10} train_model_by_chr.slurm'

geno_prefixes = list(GENOTYPE_INTER_PREFIX)
gene_annot = INTER_DIR + GENE_ANNOT_INTER2
for i, study in enumerate(STUDY_NAMES):
    expression_RDS = INTER_DIR + EXPRESSION_INTER_DIR + study + ".RDS"
    geno = INTER_DIR + GENOTYPE_INTER_DIR + geno_prefixes[0] + '.chr'
    snp_annot = INTER_DIR + SNP_ANN_INTER_DIR + SNP_ANN_INTER_PREFIX2
    cmd = CMD.format(study,expression_RDS,geno,gene_annot,snp_annot,
        N_K_FOLDS,ALPHA,MODEL_BY_CHR_DIR,SNPSET,WINDOW,HOME_DIR)
    print(cmd)
    subprocess.call(cmd, shell=True)
    time.sleep(2)
