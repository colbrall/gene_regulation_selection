#!/bin/bash

module load python/2

tiss=$1

/project/mathilab/colbranl/gene_regulation_selection/bin/PrediXcan.py --predict --dosages data/1kG_full_dosage/ --dosages_prefix chr --samples sample.txt --weights data/zhou2020_JTI_models/${tiss}.db --output_dir results/1kG_v8_JTI_predictions/${tiss}/
