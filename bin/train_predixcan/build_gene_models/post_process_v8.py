#!/usr/bin/env python

# updated to run on Slurm with Python 3 by Laura Colbran 12/2018
# first don't forget to run:
# ml GCC/5.4.0-2.26  OpenMPI/1.10.3 R/3.3.3 R-bundle-Bioconductor

import subprocess

from model_parameters_v8_1kGvariants import *
# from model_parameters_top500k import *
# from model_parameters_v8_all import *
# from model_parameters_1240k import *

for i, study in enumerate(STUDY_NAMES):
    subprocess.call(
        ['../data_processing/make_all_results_v8.sh',study,ALL_RESULTS_FILES[i],ALPHA,SNPSET])
    subprocess.call(
        ['../data_processing/make_all_betas_v8.sh',study,ALL_BETAS_FILES[i],ALPHA,SNPSET])
    subprocess.call(['../data_processing/make_all_logs_v8.sh',study,ALL_LOGS_FILES[i]])
    subprocess.call(
        ['../data_processing/make_all_covariances_v8.sh',study,ALL_COVARIANCES_FILES[i],ALPHA,SNPSET])

for i, study, in enumerate(STUDY_NAMES):
    cmd = '../data_processing/make_sqlite_db.py --output {0} --betas {1} --results {2} --construction {3} --meta {4}'.format(
        DB_FILES[i], ALL_BETAS_FILES[i], ALL_RESULTS_FILES[i], ALL_LOGS_FILES[i], ALL_META_DATA_FILES[i])
    subprocess.call(cmd, shell=True)

for i, study in enumerate(STUDY_NAMES):
    print("Filtering " + study + " on significance.")
    subprocess.call(['Rscript', '../data_processing/filter_on_significance.R', DB_FILES[i],
        INTER_DIR + GENE_ANNOT_INTER2, FILTERED_DB_FILES[i]])
