#!/bin/sh

module load R3env/3.3/3.3

export RUN=171214_D00390_0336_ACBG26ANXX
export GBS_BIN=/dataset/hiseq/active/bin/DECONVQC
export BUILD_ROOT=/dataset/hiseq/scratch/postprocessing

#psql -U agrbrdf -d agrbrdf -h invincible -f $GBS_BIN/database/extract_sample_species.psql
#$GBS_BIN/summarise_global_hiseq_taxonomy.sh $RUN
#./summarise_global_hiseq_taxonomy.sh $RUN
Rscript --vanilla ./taxonomy_clustering.r run_name=$RUN
cp -p /dataset/hiseq/scratch/postprocessing/euk_taxonomy_clustering_${RUN}.jpg ./figure_5.jpg
cp -p /dataset/hiseq/scratch/postprocessing/Clustering-of-*${RUN}*.txt .

