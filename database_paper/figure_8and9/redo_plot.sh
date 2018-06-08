#!/bin/sh

function draw_cattle() {
# note that you need to edit the R code to comment out / in  cattle/ ryegrass labels
RUN=160531_D00390_0253_AC9AB9ANXX
/dataset/bioinformatics_dev/active/R3.3/R-3.3.0/bin/Rscript --vanilla  kmer_plots_gbs.r datafolder=/dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ0223.processed_sample/uneak/kmer_analysis 1>plots.stdout 2>plots.stderr
cp -p /dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ0223.processed_sample/uneak/kmer_analysis/zipfian_distances.jpg ./Figure9.jpg
cp -p /dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ0223.processed_sample/uneak/kmer_analysis/kmer_zipfian_comparisons.jpg ./Figure8_left.jpg
}


function draw_ryegrass() {
# note that you need to edit the R code to comment out / in  cattle/ ryegrass labels
RUN=160810_D00390_0261_BC9NAVANXX
/dataset/bioinformatics_dev/active/R3.3/R-3.3.0/bin/Rscript --vanilla  kmer_plots_gbs.r datafolder=/dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ2575.processed_sample/uneak/kmer_analysis 1>plots.stdout 2>plots.stderr
cp -p /dataset/hiseq/scratch/postprocessing/$RUN.gbs/SQ2575.processed_sample/uneak/kmer_analysis/kmer_zipfian_comparisons.jpg ./Figure8_right.jpg
}

#draw_cattle
draw_ryegrass

