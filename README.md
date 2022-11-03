# gene_regulation_selection

code for Qx-related analyses on PrediXcan models

MANIFEST:
1kG_expression.jl
  calculates summary stats and makes plots for observed vs predicted expression in 1kG. Fig 1A and 1B
af_heatmap.jl
  plots allele frequency heatmaps. eg Fig 4
ancestral_match.R
  Lin's script to polarize variants for iHS
best_models.jl
  identifies the model with highest R2 for each gene
compare_rank.jl
  calculates rank correlations for p-values vs both technical and selection metrics. Fig 3A
fdr_gene_heatmap.jl
  plots the median pred/obs expression heatmaps for the top genes. eg Fig 3B
filter_dbs.jl
  given output form best_models.jl, builds a db with those specific models for convenience.
interpolate_map.R
  Iain's script to interpolate recombination maps for iHS
join_selection_output.jl
  combines selscan output into tab-delim file for downstream analyses
match_snps.jl
  for each JTI model SNP, selects a random number in the same AF bin in GTEx
PrediXcan.py
  given dosage files and a db, predicts gene expression for each individual
predixcan_search.jl
  searches dbs for genes or SNPs
predixcan.sh
  bash script to run PrediXcan.py
qqplot.jl
  plots qqplots. Fig 2
qx_stats.jl
  --empirical_p flag used to calculate gamma-corrected p-values
qx_with_p.jl
  calculates Qx statistic, and also calculates permutation-based p-values
selscan_ihs.sh
  bash script to calculate iHS using selscan (https://github.com/szpiech/selscan)
selscan.sh
  bash script to calculate nSL using selscan (https://github.com/szpiech/selscan)
sel_summary.jl
  summarizes selection stats for each gene
stitch_afs.jl
  stitches VCFtools frequency output files into combined files for use in Qx calculation
tiered_enrichment.jl
  calculates enrichment and p-values across a series of thresholds
transpose_pred_output.jl
  transposes output of PrediXcan.py
vcf2dosage.py
