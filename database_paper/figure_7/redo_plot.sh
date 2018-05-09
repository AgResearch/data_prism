#!/bin/sh

/dataset/bioinformatics_dev/active/R3.3/R-3.3.0/bin/Rscript --vanilla  kmer_entropy.r
mv kmer_zipfian_comparisons.jpg figure_7.jpg

